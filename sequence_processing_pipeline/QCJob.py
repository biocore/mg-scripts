from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os import walk, stat, listdir, makedirs
from os.path import exists, join, split
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from shutil import move
import logging
import re

logging.basicConfig(level=logging.DEBUG)


class QCJob(Job):
    def __init__(self, fastq_root_dir, output_path, sample_sheet_path,
                 minimap_database_paths, kraken2_database_path, queue_name,
                 node_count, nprocs, wall_time_limit, jmem, fastp_path,
                 minimap2_path, samtools_path, modules_to_load, qiita_job_id,
                 pool_size, max_array_length):
        """
        Submit a Torque job where the contents of fastq_root_dir are processed
        using fastp, minimap, and samtools. Human-genome sequences will be
        filtered out, if needed.
        :param fastq_root_dir: Path to a dir of Fastq files, org. by project.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: Path to a sample sheet file.
        :param minimap_database_paths: Path to human genome databases in env.
        :param kraken2_database_path: Path to human genome dat abase in env.
        :param queue_name: Torque queue name to use in running env.
        :param node_count: Number of nodes to use in running env.
        :param nprocs: Number of processes to use in runing env.
        :param wall_time_limit: Hard wall-clock-time limit for processes.
        :param jmem: String representing total memory limit for entire job.
        :param fastp_path: The path to the fastp executable
        :param minimap2_path: The path to the minimap2 executable
        :param samtools_path: The path to the samtools executable
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        :param pool_size: The number of jobs to process concurrently.
        """
        super().__init__(fastq_root_dir,
                         output_path,
                         'QCJob',
                         [fastp_path, minimap2_path, samtools_path],
                         max_array_length,
                         modules_to_load=modules_to_load)
        self.sample_sheet_path = sample_sheet_path
        self._file_check(self.sample_sheet_path)
        metadata = self._process_sample_sheet()
        self.sample_ids = metadata['sample_ids']
        self.project_data = metadata['projects']
        self.needs_trimming = metadata['needs_adapter_trimming']
        self.nprocs = 16 if nprocs > 16 else nprocs
        self.chemistry = metadata['chemistry']
        self.minimap_database_paths = minimap_database_paths
        self.kraken2_database_path = kraken2_database_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.jmem = jmem
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size
        self.suffix = 'fastq.gz'
        self.counts = {}

        self.minimum_bytes = 3100

        self.script_paths = {}

        if not isinstance(self.needs_trimming, bool):
            raise ValueError("needs_adapter_trimming must be boolean.")

        for project in self.project_data:
            project_name = project['Sample_Project']
            fastq_files = self._find_fastq_files(project_name)
            if not fastq_files:
                raise PipelineError("Fastq files could not be found for "
                                    f"project '{project_name}'")

            if project['ForwardAdapter'] == 'NA':
                project['ForwardAdapter'] = None

            if project['ReverseAdapter'] == 'NA':
                project['ReverseAdapter'] = None

            if project['ForwardAdapter'] is None:
                if project['ReverseAdapter'] is not None:
                    raise ValueError(("ForwardAdapter is declared but not "
                                      "ReverseAdapter."))

            if project['ReverseAdapter'] is None:
                if project['ForwardAdapter'] is not None:
                    raise ValueError(("ReverseAdapter is declared but not "
                                      "ForwardAdapter."))

            if not isinstance(project['HumanFiltering'], bool):
                raise ValueError("needs_adapter_trimming must be boolean.")

            spath, cnt = self._generate_job_script(project_name,
                                                   project['ForwardAdapter'],
                                                   project['ReverseAdapter'],
                                                   project['HumanFiltering'],
                                                   fastq_files)

            self.script_paths[project_name] = spath
            self.counts[project_name] = cnt

    def _filter(self, filtered_directory, empty_files_directory,
                minimum_bytes):
        empty_list = []

        for entry in listdir(filtered_directory):
            if '_R1_' in entry:
                reverse_entry = entry.replace('_R1_', '_R2_')
                full_path = join(filtered_directory, entry)
                full_path_reverse = join(filtered_directory, reverse_entry)
                if stat(full_path).st_size <= minimum_bytes or stat(
                        full_path_reverse).st_size <= minimum_bytes:
                    logging.debug(f'moving {entry} and {reverse_entry}'
                                  f' to empty list.')
                    empty_list.append(full_path)
                    empty_list.append(full_path_reverse)

        if empty_list:
            logging.debug(f'making directory {empty_files_directory}')
            makedirs(empty_files_directory, exist_ok=True)

        for item in empty_list:
            logging.debug(f'moving {item}')
            move(item, empty_files_directory)

    def _was_successful(self, project_name):
        completed_files = self._find_files(self.output_path)
        completed_files = [x for x in completed_files if
                           x.endswith('.completed') and project_name in x]

        if len(completed_files) == self.counts[project_name]:
            return True

        return False

    def run(self, callback=None):
        for project in self.project_data:
            project_name = project['Sample_Project']
            needs_human_filtering = project['HumanFiltering']
            pbs_job_id = self.qsub(self.script_paths[project_name], None, None,
                                   exec_from=self.log_path, callback=callback)
            logging.debug(f'QCJob {pbs_job_id} completed')

            source_dir = join(self.output_path, project_name)

            if not self._was_successful(project_name):
                raise PipelineError("QCJob did not complete successfully.")

            if needs_human_filtering is True:
                filtered_directory = join(source_dir, 'filtered_sequences')
            else:
                filtered_directory = join(source_dir, 'trimmed_sequences')
            empty_files_directory = join(source_dir, 'zero_files')
            self._filter(filtered_directory, empty_files_directory,
                         self.minimum_bytes)

    def _process_sample_sheet(self):
        sheet = KLSampleSheet(self.sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % self.sample_sheet_path
            raise PipelineError(s)

        header = valid_sheet.Header
        chemistry = header['chemistry']
        needs_adapter_trimming = (True if
                                  header['Assay'] == 'Metagenomics'
                                  else False)

        sample_ids = []
        for sample in valid_sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        bioinformatics = valid_sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.
        lst = bioinformatics.to_dict('records')

        # convert true/false and yes/no strings to true boolean values.
        for record in lst:
            for key in record:
                if record[key].strip().lower() in ['true', 'yes']:
                    record[key] = True
                elif record[key].strip().lower() in ['false', 'no']:
                    record[key] = False

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry,
                'projects': lst,
                'sample_ids': sample_ids,
                'needs_adapter_trimming': needs_adapter_trimming
                }

    def _find_fastq_files_in_root_dir(self, project_name):
        # project_names defined in sample-sheets must now end in a qiita id.
        # many legacy run-directories will include sub-directories named for
        # projects without this requirement. To process legacy run
        # directories, we'll include in the search-path project_names w/out
        # an id.

        search_paths = [join(self.root_dir, project_name),
                        join(self.root_dir,
                             re.sub(r'_\d+', r'', project_name)),
                        join(self.root_dir, 'Data', 'Fastq', project_name),
                        join(self.root_dir, 'Data', 'Fastq',
                             re.sub(r'_\d+', r'', project_name))]

        logging.debug("SEARCH PATHS: %s" % search_paths)

        lst = []
        for search_path in search_paths:
            if exists(search_path):
                for root, dirs, files in walk(search_path):
                    for some_file in files:
                        some_path = join(search_path, some_file)
                        if some_file.endswith('fastq.gz'):
                            lst.append(some_path)
                break
            else:
                logging.debug("PATH DOES NOT EXIST: %s" % search_path)

        # caller expects an empty list if no files were found.
        return lst

    def _find_fastq_files(self, project_name):
        # filter the list of (sample_id, sample_project) tuples stored in
        # self.sample_ids so that only the ids matching project_name are in
        # the list.
        sample_ids = filter(lambda c: c[1] == project_name, self.sample_ids)
        # strip out the project name from the matching elements.
        sample_ids = [x[0] for x in sample_ids]

        # Sample-sheet contains sample IDs, but not actual filenames.
        # Generate a list of possible fastq files to process.
        # Filter out ones that don't contain samples mentioned in sample-sheet.
        files_found = self._find_fastq_files_in_root_dir(project_name)

        lst = []

        for some_file in files_found:
            file_path, file_name = split(some_file)
            # for now, include a file if any part of its name matches any
            # sample_id w/the idea that including too many files is better
            # than silently omitting one.
            for sample_id in sample_ids:
                if sample_id in file_name:
                    lst.append(some_file)
                    break

        # caller expects an empty list if no files were found.
        return lst

    def _generate_job_script(self, project_name, adapter_a, adapter_A,
                             h_filter, fastq_file_paths):
        lines = []

        details_file_name = f'{self.job_name}_{project_name}.array-details'
        sh_details_fp = join(self.output_path, details_file_name)

        cmds = self._generate_commands(fastq_file_paths, project_name,
                                       join(self.output_path, project_name),
                                       h_filter, adapter_a, adapter_A)
        cmds = self._group_commands(cmds)

        step_count = len(cmds)

        lines.append("#!/bin/bash")
        job_name = f'{self.qiita_job_id}_{self.job_name}_{project_name}'

        lines.append(f"#PBS -N {job_name}")
        lines.append("#PBS -q %s" % self.queue_name)
        lines.append("#PBS -l nodes=%d:ppn=%d" % (self.node_count,
                                                  self.nprocs))
        lines.append("#PBS -V")
        lines.append("#PBS -l walltime=%d:00:00" % self.wall_time_limit)
        lines.append(f"#PBS -l mem={self.jmem}")

        lines.append("#PBS -t 1-%d%%%d" % (step_count, self.pool_size))

        lines.append("set -x")
        lines.append('date')
        lines.append('hostname')
        lines.append('echo ${PBS_JOBID} ${PBS_ARRAYID}')
        lines.append(f'cd {self.output_path}')

        if self.modules_to_load:
            lines.append("module load " + ' '.join(self.modules_to_load))

        lines.append('offset=${PBS_ARRAYID}')
        lines.append('step=$(( $offset - 0 ))')
        lines.append(f'cmd0=$(head -n $step {sh_details_fp} | tail -n 1)')
        lines.append('eval $cmd0')

        sentinel_file = f'{self.job_name}_{project_name}_$step.completed'
        lines.append(f'echo "Cmd Completed: $cmd0\n" > logs/{sentinel_file}')

        job_script_path = join(self.output_path,
                               f'{self.job_name}_{project_name}.sh')

        with open(job_script_path, 'w') as f:
            f.write('\n'.join(lines))

        with open(sh_details_fp, 'w') as f:
            f.write('\n'.join(cmds))

        return job_script_path, step_count

    def _generate_commands(self, fastq_file_paths, project_name, project_dir,
                           h_filter, adapter_a, adapter_A):
        """
        Generate the command-lines needed to QC the data, based on the
        parameters supplied to the object. The result will be the list of
        strings needed to process the fastq files referenced in the trim file.
        :return: A list of strings that can be run using Popen() and the like.
        """
        fastp_reports_dir = join(project_dir, 'fastp_reports_dir')

        # if this directory is not made, then fastp will not create the html
        # and json directories and generate output files for them.
        makedirs(fastp_reports_dir, exist_ok=True)

        cmds = []

        od = 'filtered_sequences' if h_filter is True else 'trimmed_sequences'

        if self.needs_trimming is True:
            if not exists(join(project_dir, od)):
                makedirs(join(project_dir, od),
                         exist_ok=True)

        if self.needs_trimming is True:
            paths = [x.strip() for x in fastq_file_paths if '_R1_' in x]
            for fastq_file_path in paths:
                current_dir = split(fastq_file_path)[0]
                _, filename1 = split(fastq_file_path)
                filename2 = filename1.replace('_R1_00', '_R2_00')

                if h_filter is True:
                    cmds.append(self._gen_chained_cmd(current_dir,
                                                      filename1,
                                                      filename2,
                                                      fastp_reports_dir,
                                                      project_name,
                                                      project_dir,
                                                      adapter_a,
                                                      adapter_A))
                else:
                    cmds.append(self._gen_fastp_cmd(current_dir,
                                                    filename1,
                                                    filename2,
                                                    project_dir,
                                                    adapter_a,
                                                    adapter_A))
            return cmds

        logging.warning("QCJob created w/a_trim set to False.")
        return None

    def _gen_fastp_cmd(self, current_dir, filename1, filename2, project_dir,
                       adapter_a, adapter_A):
        """
        Generates a command-line string for running fastp, based on the
         parameters supplied to the object.
        :return: A string suitable for executing in Popen() and the like.
        """
        read1_input_path = join(current_dir, filename1)
        read2_input_path = join(current_dir, filename2)

        partial_path = join(current_dir, project_dir)
        json_output_path = join(partial_path, 'json',
                                filename1.replace('.fastq.gz', '.json'))
        html_output_path = join(partial_path, 'html',
                                filename1.replace('.fastq.gz', '.html'))
        read1_output_path = join(partial_path, 'trimmed_sequences',
                                 filename1.replace('.fastq.gz',
                                                   '.fastp.fastq.gz'))
        read2_output_path = join(partial_path, 'trimmed_sequences',
                                 filename2.replace('.fastq.gz',
                                                   '.fastp.fastq.gz'))
        report_title = filename1.replace('.fastq.gz', '') + '_report'

        result = self.fastp_path
        if adapter_a:
            # assume that if adapter_a is not None, then adapter_A is not None
            # as well. We performed this check already.
            result += (f' --adapter_sequence {adapter_a}'
                       f' --adapter_sequence_r2 {adapter_A}')

        result += (f' -l 100 -i {read1_input_path} -I {read2_input_path} -w '
                   f'{self.nprocs} -j {json_output_path} -h {html_output_path}'
                   f' -o {read1_output_path} -O {read2_output_path} -R '
                   f'{report_title}')

        return result

    def _gen_chained_cmd(self, current_dir, filename1, filename2,
                         fastp_reports_dir, project_name, project_dir,
                         adapter_a, adapter_A):
        read1_input_path = join(current_dir, filename1)
        read2_input_path = join(current_dir, filename2)

        tmp_path = join(project_name, fastp_reports_dir, 'json')
        makedirs(tmp_path, exist_ok=True)
        json_output_path = join(tmp_path, filename1.replace('.fastq.gz',
                                                            '.json'))

        tmp_path = join(project_name, fastp_reports_dir, 'html')
        makedirs(tmp_path, exist_ok=True)
        html_output_path = join(tmp_path, filename1.replace('.fastq.gz',
                                                            '.html'))

        partial = join(project_dir, 'filtered_sequences')

        path1 = join(partial, filename1.replace('.fastq.gz',
                                                '.trimmed.fastq.gz'))
        path2 = join(partial, filename2.replace('.fastq.gz',
                                                '.trimmed.fastq.gz'))

        kraken_report_path = join(partial,
                                  filename1.replace('.fastq.gz',
                                                    '.kraken2_report.txt'))

        kraken_output_path = join(partial,
                                  filename1.replace('.fastq.gz',
                                                    '.kraken2.trimmed.#.fastq')
                                  )

        result = self.fastp_path

        if adapter_a:
            # assume that if adapter_a is not None, then adapter_A is not None
            # as well. We performed this check already.
            result += (f' --adapter_sequence {adapter_a}'
                       f' --adapter_sequence_r2 {adapter_A}')

        result += (f' -l 100 -i {read1_input_path} -I {read2_input_path} -w '
                   f'{self.nprocs} -j {json_output_path} -h {html_output_path}'
                   ' --stdout ')

        # for each database specified in configuration, run the data through
        # minimap using the database and pipe to samtools to tame the output.
        for database in self.minimap_database_paths:
            result += (f'| {self.minimap2_path} -ax sr -t {self.nprocs} '
                       f'{database} - -a | {self.samtools_path} fastq -@ '
                       f'{self.nprocs} -f 12 -F 256 ')

        # append the final parameters to write out the final output to disk
        result += f'-1 {path1} -2 {path2}'

        # add krakken2
        result += (f'\nkraken2 --threads {self.nprocs} --db '
                   f'{self.kraken2_database_path} --report '
                   f'{kraken_report_path} --unclassified-out '
                   f'{kraken_output_path} --paired {path1} {path2}')

        return result
