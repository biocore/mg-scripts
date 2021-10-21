from os import listdir, makedirs
from os.path import join, basename, exists
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from functools import partial
import logging


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

        self.project_names = []
        self.commands = self._get_commands()
        split_count = self._generate_split_count(len(self.commands))
        lines_per_split = int((len(self.commands) + split_count - 1) /
                              split_count)
        self._generate_trim_files(self.commands, lines_per_split)
        self.script_path = self._generate_job_script()

        for project_name in self.project_names:
            p_path = partial(join, self.products_dir, 'FastQC', project_name)
            makedirs(p_path('bclconvert'), exist_ok=True)
            makedirs(p_path('filtered_sequences'), exist_ok=True)

    def _get_commands(self):
        results = []

        # gather the parameters for processing all relevant raw fastq files.
        params = self._scan_fastq_files(True)

        for thread_count, fwd_file_path, rev_file_path, output_path in params:
            command = ['fastqc', '--noextract', '-t', str(thread_count),
                       fwd_file_path, rev_file_path, '-o', output_path]
            results.append(' '.join(command))

        # next, do the same for the trimmed/filtered fastq files.
        params = self._scan_fastq_files(False)

        for thread_count, fwd_file_path, rev_file_path, output_path in params:
            command = ['fastqc', '--noextract', '-t', str(thread_count),
                       fwd_file_path, rev_file_path, '-o', output_path]
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

            # break files up into R1, R2, I1, I2
            # assume _R1_ does not occur in the path as well.
            r1_only = [x for x in files if '_R1_' in x]
            r2_only = [x for x in files if '_R2_' in x]

            if len(r1_only) != len(r2_only):
                raise PipelineError('counts of R1 and R2 files do not match')

            i1_only = [x for x in files if '_I1_' in x]
            i2_only = [x for x in files if '_I2_' in x]

            if len(i1_only) != len(i2_only):
                raise PipelineError('counts of I1 and I2 files do not match')

            if r1_only:
                tmp = ' '.join(r1_only)
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

                r1_only.sort()
                r2_only.sort()

                if filter_type == 'amplicon':
                    i1_only.sort()
                    i2_only.sort()
                    results.append((directory, filter_type, project_dir,
                                    r1_only + i1_only, r2_only + i2_only))

                results.append(
                    (directory, filter_type, project_dir, r1_only, r2_only))

        return results

    def _scan_fastq_files(self, is_raw_input=False):
        find_path = (self.raw_fastq_path if is_raw_input else
                     self.processed_fastq_path)

        projects = self._find_projects(find_path, is_raw_input)

        fastqc_results = []

        for proj_name, fltr_type, fastq_fp, fwd_files, rev_files in projects:
            self.project_names.append(proj_name)
            base_path = join(self.fastqc_output_path, proj_name)
            dir_name = 'bclconvert' if is_raw_input else fltr_type

            for some_fwd_file, some_rev_file in zip(fwd_files, rev_files):
                fastqc_results.append((self.nthreads, some_fwd_file,
                                       some_rev_file,
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
            # MultiQC doesn't like input paths that don't exist. Simply add
            # all paths that do exist as input.
            input_path_list = []
            p_path = partial(join, self.run_dir, self.run_id, 'FastQC',
                             project)

            for filter_type in ['bclconvert', 'trimmed_sequences',
                                'filtered_sequences', 'amplicon']:
                input_path_list.append(p_path(filter_type))

            input_path_list.append(join(self.run_dir, self.run_id, 'FastQC',
                                        'Reports'))

            input_path_list.append(join(self.processed_fastq_path, project,
                                        'fastp_reports_dir', 'json'))

            # I don't usually see a json directory associated with raw data.
            # It looks to be metadata coming directly off the machine, in the
            # few instances I've seen it in /sequencing...
            input_path_list.append(join(self.raw_fastq_path, project, 'json'))

            input_path_list = [x for x in input_path_list if exists(x)]

            cmd_head = ['multiqc', '-c', self.multiqc_config_file_path,
                        '--fullnames', '--force']

            cmd_tail = ['-o', join(self.fastqc_output_path, project,
                                   'multiqc'), '--interactive']

            results = self._system_call(' '.join(cmd_head + input_path_list
                                                 + cmd_tail))

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
        if self.modules_to_load:
            lines.append("module load " + ' '.join(self.modules_to_load))
        lines.append('offset=${PBS_ARRAYID}')
        lines.append('step=$(( $offset - 0 ))')
        lines.append(f'cmd0=$(head -n $step {sh_details_fp} | tail -n 1)')
        lines.append('eval $cmd0')

        with open(job_script_path, 'w') as f:
            f.write('\n'.join(lines))

        with open(sh_details_fp, 'w') as f:
            f.write('\n'.join(self.commands))

        return job_script_path
