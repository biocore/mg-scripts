from sequence_processing_pipeline.Job import Job
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join, split, basename
from os import walk, remove, stat, listdir, makedirs
import logging
from sequence_processing_pipeline.QCHelper import QCHelper
from shutil import move


class QCJob(Job):
    def __init__(self, run_dir, sample_sheet_path, mmi_db_path, queue_name,
                 node_count, nprocs, wall_time_limit, jmem, fastp_path,
                 minimap2_path, samtools_path, modules_to_load, qiita_job_id,
                 pool_size, products_dir):
        '''
        Submit a Torque job where the contents of run_dir are processed using
        fastp, minimap, and samtools. Human-genome sequences will be filtered
        out if needed.
        :param run_dir: Path to a run directory.
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
        :param products_dir: The path to the products directory
        '''
        # for now, keep this run_dir instead of abspath(run_dir)
        self.job_name = 'QCJob'
        super().__init__(run_dir,
                         self.job_name,
                         [fastp_path, minimap2_path, samtools_path],
                         modules_to_load)
        self.sample_sheet_path = sample_sheet_path
        metadata = self._process_sample_sheet()
        self.sample_ids = metadata['sample_ids']
        self.project_data = metadata['projects']
        self.needs_a_trimming = metadata['needs_adapter_trimming']
        self.trim_file = 'split_file_'
        self.nprocs = nprocs
        self.chemistry = metadata['chemistry']
        self.mmi_db_path = mmi_db_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.nprocs = nprocs
        self.wall_time_limit = wall_time_limit
        self.run_id = basename(self.run_dir)
        self.products_dir = products_dir
        self.jmem = jmem
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        self.modules_to_load = modules_to_load
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size

        # set to 500 bytes to avoid empty and small files that Qiita
        # has trouble with.
        self.minimum_bytes = 500

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

    def run(self):
        for project in self.project_data:
            fastq_files = self._find_fastq_files(project['Sample_Project'])
            split_count = self._generate_split_count(len(fastq_files))
            self._clear_trim_files()
            lines_per_split = int((
                len(fastq_files) + split_count - 1) / split_count)
            trim_files = self._generate_trim_files(
                fastq_files, lines_per_split)
            if len(trim_files) != split_count:
                logging.warning("The number of trim files does not equal the "
                                "number of expected files.")

            script_path = self._generate_job_script(project['Sample_Project'],
                                                    project['ForwardAdapter'],
                                                    project['ReverseAdapter'],
                                                    project['HumanFiltering'])

            pbs_job_id = self.qsub(script_path, None, None)
            logging.debug(f'QCJob {pbs_job_id} completed')

            source_dir = join(self.products_dir, project['Sample_Project'])
            filtered_directory = join(source_dir, 'filtered_sequences')
            empty_files_directory = join(source_dir, 'zero_files')
            self._filter(filtered_directory, empty_files_directory,
                         self.minimum_bytes)

    def _generate_trim_files(self, fastq_files, split_count):
        def _chunk_list(some_list, n):
            # taken from https://bit.ly/38S9O8Z
            for i in range(0, len(some_list), n):
                yield some_list[i:i + n]

        new_files = []
        count = 0
        for chunk in _chunk_list(fastq_files, split_count):
            trim_file_name = '%s%d' % (self.trim_file, count)
            trim_file_path = join(self.run_dir, trim_file_name)
            with open(trim_file_path, 'w') as f:
                for line in chunk:
                    f.write("%s\n" % line)
            new_files.append(trim_file_path)
            count += 1

        return new_files

    def _clear_trim_files(self):
        # remove all files with a name beginning in self.trim_file.
        # assume cleaning the entire run_dir is overkill, but won't
        # hurt anything.
        for root, dirs, files in walk(self.run_dir):
            for some_file in files:
                if self.trim_file in some_file:
                    some_path = join(root, some_file)
                    remove(some_path)

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
            for key in ['PolyGTrimming', 'HumanFiltering']:
                if record[key].strip().lower() in ['true', 'yes']:
                    record[key] = True
                elif record[key].strip().lower() in ['false', 'no']:
                    record[key] = False
                else:
                    raise PipelineError(f"Unexpected value '{record[key]}' "
                                        f"for column '{key}' in sample-sheet '"
                                        f"{self.sample_sheet_path}'")

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry,
                'projects': lst,
                'sample_ids': sample_ids,
                'needs_adapter_trimming': needs_adapter_trimming
                }

    def _find_fastq_files_in_run_dir(self, project_name):
        search_path = join(self.run_dir, 'Data', 'Fastq', project_name)
        lst = []
        for root, dirs, files in walk(search_path):
            for some_file in files:
                if some_file.endswith('fastq.gz'):
                    some_path = join(search_path, some_file)
                    lst.append(some_path)
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
        files_found = self._find_fastq_files_in_run_dir(project_name)
        logging.debug(f"Fastq files found: {files_found}")

        lst = []

        for some_file in files_found:
            file_path, file_name = split(some_file)
            # assume each file_name begins with a sample_id.
            # This name is generated by programs like bcl-convert and
            # bcl2fastq so this is fairly durable.
            # e.g.: X00185910 and X00185910_S177_L003_R1_001.fastq.gz.
            sample_id = file_name.split('_')[0]
            if sample_id in sample_ids:
                # this Fastq file is one we should process.
                logging.debug(f"Located Fastq file for id '{id}': {some_file}")
                # Save the full path.
                lst.append(some_file)

        return lst

    def _generate_split_count(self, count):
        if count > 2000:
            return 16
        elif count <= 2000 and count > 1000:
            return 10
        elif count <= 1000 and count > 500:
            return 4

        return 1

    def _generate_job_script(self, project_name, adapter_a, adapter_A,
                             h_filter):
        lines = []

        # Unlike w/ConvertBCL2FastqJob, multiple job scripts are generated,
        # one for each project defined in the sample sheet.
        # This Job() class does not use Job.stdout_log_path and
        # Job.stderr_log_path. Instead, it uses numbered paths provided by
        # generate_job_script_path().
        job_script_path, output_log_path, error_log_path = \
            self.generate_job_script_path()

        project_products_dir = join(self.products_dir, project_name)

        tf_fp = join(self.run_dir, self.trim_file + '0')
        sh_details_fp = join(
            self.run_dir, self.trim_file + project_name + '.array-details')

        qc = QCHelper(self.nprocs, tf_fp, project_name, project_products_dir,
                      self.mmi_db_path, adapter_a, adapter_A,
                      self.needs_a_trimming, h_filter, self.chemistry,
                      self.fastp_path, self.minimap2_path, self.samtools_path)

        cmds = qc.generate_commands()

        lines.append("#!/bin/bash")
        # declare a name for this job
        lines.append("#PBS -N %s" %
                     f"{self.qiita_job_id}_QCJob_{project_name}")

        # what torque calls a queue, slurm calls a partition
        lines.append("#PBS -q %s" % self.queue_name)

        # request a single node, and n processes for this job.
        lines.append("#PBS -l nodes=%d:ppn=%d" % (self.node_count,
                                                  self.nprocs))
        lines.append("#PBS -V")
        lines.append("#PBS -l walltime=%d:00:00" % self.wall_time_limit)
        lines.append(f"#PBS -l mem={self.jmem}")

        # Generate output files
        lines.append("#PBS -o localhost:%s.${PBS_ARRAYID}" % output_log_path)
        lines.append("#PBS -e localhost:%s.${PBS_ARRAYID}" % error_log_path)

        # Configure array size, and the size of the pool of concurrent jobs.
        lines.append("#PBS -t 1-%d%%%d" % (len(cmds), self.pool_size))

        lines.append("set -x")
        lines.append('date')
        lines.append('hostname')
        lines.append('echo ${PBS_JOBID} ${PBS_ARRAYID}')
        lines.append("cd %s" % self.run_dir)
        modules_to_load = "module load " + ' '.join(self.modules_to_load)
        logging.debug(f"QCJob Modules to load: {modules_to_load}")
        lines.append(modules_to_load)
        lines.append('offset=${PBS_ARRAYID}')
        lines.append('step=$(( $offset - 0 ))')
        lines.append(f'cmd0=$(head -n $step {sh_details_fp} | tail -n 1)')
        lines.append('eval $cmd0')

        with open(job_script_path, 'w') as f:
            f.write('\n'.join(lines))

        with open(sh_details_fp, 'w') as f:
            f.write('\n'.join(cmds))

        return job_script_path
