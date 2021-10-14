from os import listdir
from os.path import join, basename
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
import logging


logging.basicConfig(level=logging.DEBUG)


class FastQCJob(Job):
    def __init__(self, run_dir, output_directory, nprocs, nthreads,
                 fastqc_path, modules_to_load, qiita_job_id, queue_name,
                 node_count, wall_time_limit, jmem, pool_size,
                 multiqc_config_file_path):
        self.job_name = 'FastQCJob'
        super().__init__(run_dir, self.job_name, [fastqc_path],
                         modules_to_load)
        self.run_id = basename(self.run_dir)
        self.nprocs = nprocs
        self.nthreads = nthreads
        self.fastqc_path = fastqc_path

        if self._file_check(multiqc_config_file_path):
            self.multiqc_config_file_path = multiqc_config_file_path
        else:
            raise PipelineError(f"'{multiqc_config_file_path}' does not exist")

        self.modules_to_load = modules_to_load
        self.raw_fastq_path = join(self.run_dir, 'Data', 'Fastq')
        self.fastqc_output_path = output_directory
        self.processed_fastq_path = join(self.run_dir, self.run_id)
        self.qiita_job_id = qiita_job_id
        self.trim_file = 'fqc_split_file_'
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.products_dir = join(self.run_dir, self.run_id)
        # for CI, create this directory if it doesn't exist already.
        self._directory_check(self.products_dir, create=True)

        self.jmem = jmem
        self.pool_size = pool_size

        self.project_names = []  # ['THDMI_US_10317']
        self.commands = self._get_commands()
        split_count = self._generate_split_count(len(self.commands))
        lines_per_split = int((len(self.commands) + split_count - 1) /
                              split_count)
        self._generate_trim_files(self.commands, lines_per_split)
        self.script_path = self._generate_job_script()

    def _get_commands(self):
        results = []

        # gather the parameters for processing all relevant raw fastq files.
        params = self._scan_fastq_files(True)

        for number_of_threads, file_path, output_path in params:
            command = ['fastqc', '--noextract', '-t', str(number_of_threads),
                       file_path, '-o', output_path]
            results.append(' '.join(command))

        # next, do the same for the trimmed/filtered fastq files.
        params = self._scan_fastq_files(False)

        for number_of_threads, file_path, output_path in params:
            command = ['fastqc', '--noextract', '-t', str(number_of_threads),
                       file_path, '-o', output_path]
            results.append(' '.join(command))

        # remove duplicate project names from the list, now that
        # processing is complete.
        self.project_names = list(set(self.project_names))

        return results

    def _find_projects(self, path_to_run_id_data_fastq_dir, is_raw_input):
        results = []
        for directory in listdir(path_to_run_id_data_fastq_dir):
            project_dir = join(path_to_run_id_data_fastq_dir, directory)
            files = self._find_files(project_dir)

            # extract only fastq files from the list
            files = [x for x in files if x.endswith('.fastq.gz')]

            if files:
                tmp = ' '.join(files)
                if 'trimmed_sequences' in tmp:
                    # a_trim = True, h_filter= = False
                    filter_type = 'trimmed_sequences'
                elif 'filtered_sequences' in tmp:
                    # a_trim = True, h_filter= = True
                    # (trimmed AND filtered)
                    filter_type = 'filtered_sequences'
                elif 'amplicon' in tmp:
                    # a_trim = False, h_filter= = False
                    filter_type = 'amplicon'
                else:
                    if is_raw_input:
                        filter_type = 'raw'
                    else:
                        raise ValueError("indeterminate type")

                if filter_type != 'amplicon':
                    # filter out index '_In_' files
                    files = [x for x in files if '_R1_' in x or '_R2_' in x]
                results.append((directory, filter_type, project_dir, files))

        return results

    def _scan_fastq_files(self, is_raw_input=False):
        find_path = (self.raw_fastq_path if is_raw_input else
                     self.processed_fastq_path)

        projects = self._find_projects(find_path, is_raw_input)

        fastqc_results = []

        for project_name, filter_type, fastq_path, files in projects:
            self.project_names.append(project_name)
            base_path = join(self.fastqc_output_path, project_name)
            dir_name = 'bclconvert' if is_raw_input else filter_type

            for some_file in files:
                fastqc_results.append((self.nthreads, some_file,
                                       join(base_path, dir_name)))

        return fastqc_results

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

    def _generate_split_count(self, count):
        if count > 2000:
            return 16
        elif count <= 2000 and count > 1000:
            return 10
        elif count <= 1000 and count > 500:
            return 4

        return 1

    def run(self):
        pbs_job_id = self.qsub(self.script_path, None, None)
        logging.debug(pbs_job_id)

        for project in self.project_names:
            logging.debug('PROJECT: %s' % project)
            fastp_reports = join(self.run_dir, self.run_id, project,
                                 'fastp_reports_dir', 'json')
            fastqc_reports = join(self.run_dir, self.run_id, 'FastQC',
                                  project)
            # filtered_reports may be redundant, as multiqc  can use
            # fastqc_reports as a root to search
            filtered_reports = join(fastqc_reports, 'filtered_sequences')
            cmd = ['multiqc', '-c', self.multiqc_yaml_path, '--fullnames',
                   '--force', fastp_reports, fastqc_reports, filtered_reports,
                   '-o', join(self.fastqc_output_path, 'multiqc'),
                   '--interactive']
            logging.debug(cmd)

            results = self._system_call(' '.join(cmd))
            logging.debug(f"_stdout: {results['stdout']}")
            logging.debug(f"_stderr: {results['stderr']}")
            logging.debug(f"return code: {results['return_code']}")

            if results['return_code'] != 0:
                raise PipelineError("multiqc encountered an error")

    def _generate_job_script(self):
        lines = []

        job_script_path, output_log_path, error_log_path = \
            self.generate_job_script_path()

        sh_details_fp = join(self.run_dir, self.trim_file + '.array-details')

        lines.append("#!/bin/bash")
        lines.append("#PBS -N %s" % f"{self.qiita_job_id}_FastQCJob")
        lines.append("#PBS -q %s" % self.queue_name)
        lines.append("#PBS -l nodes=%d:ppn=%d" % (self.node_count,
                                                  self.nprocs))
        lines.append("#PBS -V")
        lines.append("#PBS -l walltime=%d:00:00" % self.wall_time_limit)
        lines.append(f"#PBS -l mem={self.jmem}")
        lines.append("#PBS -o localhost:%s.${PBS_ARRAYID}" % output_log_path)
        lines.append("#PBS -e localhost:%s.${PBS_ARRAYID}" % error_log_path)
        lines.append("#PBS -t 1-%d%%%d" % (len(self.commands), self.pool_size))
        lines.append("set -x")
        lines.append('date')
        lines.append('hostname')
        lines.append('echo ${PBS_JOBID} ${PBS_ARRAYID}')
        lines.append("cd %s" % self.run_dir)
        tmp = "module load " + ' '.join(self.modules_to_load)
        logging.debug(f"QCJob Modules to load: {tmp}")
        lines.append(tmp)

        lines.append('offset=${PBS_ARRAYID}')
        lines.append('step=$(( $offset - 0 ))')
        lines.append(f'cmd0=$(head -n $step {sh_details_fp} | tail -n 1)')
        lines.append('eval $cmd0')

        with open(job_script_path, 'w') as f:
            f.write('\n'.join(lines))

        with open(sh_details_fp, 'w') as f:
            f.write('\n'.join(self.commands))

        return job_script_path
