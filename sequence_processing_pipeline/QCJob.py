from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os import walk, stat, listdir, makedirs
from os.path import exists, join, split
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.QCHelper import QCHelper
from shutil import move
import logging
import re

logging.basicConfig(level=logging.DEBUG)


class QCJob(Job):
    def __init__(self, fastq_root_dir, output_path, sample_sheet_path,
                 mmi_db_path, queue_name, node_count, nprocs, wall_time_limit,
                 jmem, fastp_path, minimap2_path, samtools_path,
                 modules_to_load, qiita_job_id, pool_size, max_array_length):
        '''
        Submit a Torque job where the contents of fastq_root_dir are processed
        using fastp, minimap, and samtools. Human-genome sequences will be
        filtered out, if needed.
        :param fastq_root_dir: Path to a dir of Fastq files, org. by project.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: Path to a sample sheet file.
        :param mmi_db_path: Path to human genome database in running env.
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
        '''
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
        self.nprocs = nprocs
        self.chemistry = metadata['chemistry']
        self.mmi_db_path = mmi_db_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.nprocs = nprocs
        self.wall_time_limit = wall_time_limit
        self.jmem = jmem
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size
        self.suffix = 'fastq.gz'

        # set to 500 bytes to avoid empty and small files that Qiita
        # has trouble with.
        self.minimum_bytes = 500

        self.script_paths = {}

        for project in self.project_data:
            project_name = project['Sample_Project']
            fastq_files = self._find_fastq_files(project_name)
            if not fastq_files:
                raise PipelineError("Fastq files could not be found for "
                                    f"project '{project_name}'")

            script_path = self._generate_job_script(project_name,
                                                    project['ForwardAdapter'],
                                                    project['ReverseAdapter'],
                                                    project['HumanFiltering'],
                                                    fastq_files)

            self.script_paths[project_name] = script_path

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

    def run(self, callback=None):
        for project in self.project_data:
            project_name = project['Sample_Project']
            needs_human_filtering = project['HumanFiltering']
            pbs_job_id = self.qsub(self.script_paths[project_name], None, None,
                                   exec_from=self.log_path, callback=callback)
            logging.debug(f'QCJob {pbs_job_id} completed')
            source_dir = join(self.output_path, project_name)
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

        project_dir = join(self.output_path, project_name)

        details_file_name = f'{self.job_name}_{project_name}.array-details'
        sh_details_fp = join(self.output_path, details_file_name)

        qc = QCHelper(self.nprocs, fastq_file_paths, project_name,
                      project_dir, self.mmi_db_path, adapter_a,
                      adapter_A, self.needs_trimming, h_filter,
                      self.chemistry, self.fastp_path, self.minimap2_path,
                      self.samtools_path)

        cmds = qc.generate_commands()
        cmds = self._group_commands(cmds)

        lines.append("#!/bin/bash")
        job_name = f'{self.qiita_job_id}_{self.job_name}_{project_name}'

        lines.append(f"#PBS -N {job_name}")
        lines.append("#PBS -q %s" % self.queue_name)
        lines.append("#PBS -l nodes=%d:ppn=%d" % (self.node_count,
                                                  self.nprocs))
        lines.append("#PBS -V")
        lines.append("#PBS -l walltime=%d:00:00" % self.wall_time_limit)
        lines.append(f"#PBS -l mem={self.jmem}")

        lines.append("#PBS -t 1-%d%%%d" % (len(cmds), self.pool_size))

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

        job_script_path = join(self.output_path,
                               f'{self.job_name}_{project_name}.sh')

        with open(job_script_path, 'w') as f:
            f.write('\n'.join(lines))

        with open(sh_details_fp, 'w') as f:
            f.write('\n'.join(cmds))

        return job_script_path
