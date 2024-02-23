from os import listdir, makedirs
from os.path import exists, join, basename
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
from functools import partial
from json import dumps
import logging
import glob


class FastQCJob(Job):
    def __init__(self, run_dir, output_path, raw_fastq_files_path,
                 processed_fastq_files_path, nprocs, nthreads, fastqc_path,
                 modules_to_load, qiita_job_id, queue_name, node_count,
                 wall_time_limit, jmem, pool_size, multiqc_config_file_path,
                 max_array_length, is_amplicon):
        super().__init__(run_dir,
                         output_path,
                         'FastQCJob',
                         [fastqc_path],
                         max_array_length,
                         modules_to_load=modules_to_load)

        self.nprocs = nprocs
        self.nthreads = nthreads
        self.fastqc_path = fastqc_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.jmem = jmem
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size
        self.raw_fastq_files_path = raw_fastq_files_path
        self.processed_fastq_files_path = processed_fastq_files_path
        self.is_amplicon = is_amplicon

        self.job_script_path = join(self.output_path, f"{self.job_name}.sh")

        self._file_check(multiqc_config_file_path)
        self.multiqc_config_file_path = multiqc_config_file_path

        self.project_names = []
        self.commands, self.project_names = self._get_commands()
        # for lists greater than n commands, chain the extra commands,
        # distributing them evenly throughout the first n commands.
        self.commands = self._group_commands(self.commands)
        self.suffix = 'fastqc.html'

        self._generate_job_script()

    def _get_commands(self):
        """
        Generate a set of commands to execute, based on the input metadata.
        :return: A list of commands to execute w/in a job script
        """
        results = []

        if self.is_amplicon:
            # skip this step for amplicon runs since raw and processed are the
            # same file.
            project_names = []
        else:
            # gather the parameters for processing all relevant raw fastq
            # files.
            params, project_names = self._scan_fastq_files(True)

            for fwd_file_path, rev_file_path, output_path in params:
                command = ['fastqc', '--noextract', '-t', str(self.nthreads),
                           fwd_file_path, rev_file_path, '-o', output_path]
                results.append(' '.join(command))

        # next, do the same for the trimmed/filtered fastq files.
        params, additional_project_names = self._scan_fastq_files(False)

        for fwd_file_path, rev_file_path, output_path in params:
            command = ['fastqc', '--noextract', '-t', str(self.nthreads),
                       fwd_file_path, rev_file_path, '-o', output_path]
            results.append(' '.join(command))

        # remove duplicate project names from the list
        project_names = list(set(project_names + additional_project_names))

        return results, project_names

    def _find_projects(self, path_to_run_id_data_fastq_dir, is_raw_input):
        results = []
        for directory in listdir(path_to_run_id_data_fastq_dir):
            project_dir = join(path_to_run_id_data_fastq_dir, directory)
            files = self._find_files(project_dir)

            # extract only fastq files from the list and remove all files
            # in 'zero_files' sub-folder.
            files = [x for x in files if x.endswith('.fastq.gz') and
                     'zero_files' not in x]

            # remove fastq files in the only-adapter-filtered
            # folder from consideration if they are present.
            files = [x for x in files if 'only-adapter-filtered' not in x]

            # break files up into R1, R2, I1, I2
            # assume _R1_ does not occur in the path as well.
            r1_only = [x for x in files if '_R1_' in x]
            r2_only = [x for x in files if '_R2_' in x]

            # amplicon runs may or may not have an i2. this is okay.
            i1_only = [x for x in files if '_I1_' in x]
            i2_only = [x for x in files if '_I2_' in x]

            if not self.is_amplicon and len(i1_only) != len(i2_only):
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

        return results if results else None

    def _scan_fastq_files(self, is_raw_input=False):
        find_path = (self.raw_fastq_files_path if is_raw_input else
                     self.processed_fastq_files_path)

        projects = self._find_projects(find_path, is_raw_input)

        if projects is None:
            raise PipelineError("There are no fastq files for FastQCJob to "
                                "process in %s." % find_path)

        fastqc_results = []

        # rather than get a list of project_names from a sample-sheet,
        # gather a list of names from the outputs of previous jobs.
        project_names = []

        for proj_name, fltr_type, fastq_fp, fwd_files, rev_files in projects:
            project_names.append(proj_name)
            p_path = partial(join, self.output_path, 'fastqc', proj_name)
            output_dir = p_path('bclconvert' if is_raw_input else fltr_type)
            makedirs(output_dir, exist_ok=True)

            for some_fwd_file, some_rev_file in zip(fwd_files, rev_files):
                fastqc_results.append((some_fwd_file, some_rev_file,
                                       output_dir))
        # remove duplicates
        project_names = list(set(project_names))

        return fastqc_results, project_names

    def _get_failed_indexes(self, job_id):
        completed_files = self._find_files(self.output_path)
        completed_files = [x for x in completed_files if
                           x.endswith('.completed')]

        completed_indexes = []

        for completed_file in completed_files:
            # remove path and .completed extension from file-name. e.g.:
            # 'project1_0', 'project1_1', ..., 'project1_n'
            completed_file = basename(completed_file).replace('.completed', '')
            # extract the line number in the .detailed file corresponding to
            # the command used for this job
            completed_indexes.append(int(completed_file.split('_')[-1]))

        # a successfully completed job array should have a list of array
        # numbers from 0 - len(self.commands).
        all_indexes = [x for x in range(1, len(self.commands) + 1)]
        failed_indexes = list(set(all_indexes) - set(completed_indexes))
        failed_indexes.sort()

        # generate log-file here instead of in run() where it can be
        # unittested more easily.
        log_fp = join(self.output_path,
                      'logs',
                      f'failed_indexes_{job_id}.json')

        if failed_indexes:
            with open(log_fp, 'w') as f:
                f.write(dumps({'job_id': job_id,
                               'failed_indexes': failed_indexes}, indent=2))

        return failed_indexes

    def run(self, callback=None):
        try:
            job_info = self.submit_job(self.job_script_path,
                                       exec_from=self.log_path,
                                       callback=callback)
        except JobFailedError as e:
            # When a job has failed, parse the logs generated by this specific
            # job to return a more descriptive message to the user.
            info = self.parse_logs()
            # prepend just the message component of the Error.
            info.insert(0, str(e))
            raise JobFailedError('\n'.join(info))

        logging.debug(job_info)

        # If project-level reports were not needed, MultiQC could simply be
        # given the path to the run-directory itself and it will discover all
        # of the relevant data files. Confirmed that for a one-project sample-
        # sheet, this produces and equivalent report.

        for project in self.project_names:
            # MultiQC doesn't like input paths that don't exist. Simply add
            # all paths that do exist as input.
            input_path_list = []
            p_path = partial(join, self.output_path, 'fastqc')

            for filter_type in ['bclconvert', 'trimmed_sequences',
                                'filtered_sequences', 'amplicon']:
                input_path_list.append(p_path(project, filter_type))

            input_path_list.append(p_path(project, 'Reports'))

            p_path = partial(join, self.processed_fastq_files_path, project)
            input_path_list.append(p_path('fastp_reports_dir', 'json'))

            # I don't usually see a json directory associated with raw data.
            # It looks to be metadata coming directly off the machine, in the
            # few instances I've seen it in /sequencing...
            p_path = partial(join, self.raw_fastq_files_path, project)
            input_path_list.append(p_path('json'))

            input_path_list = [x for x in input_path_list if exists(x)]

            cmd_head = ['multiqc', '-c', self.multiqc_config_file_path,
                        '--fullnames', '--force']

            # --interactive graphs is set to True in MultiQC configuration
            # file and hence this switch was redunant and now removed.
            cmd_tail = ['-o', join(self.output_path, 'multiqc', project)]

            cmd = ' '.join(cmd_head + input_path_list + cmd_tail)

            results = self._system_call(cmd, callback=callback)

            if results['return_code'] != 0:
                raise PipelineError("multiqc encountered an error")

        if self._get_failed_indexes(job_info['job_id']):
            # raise error if list isn't empty.
            raise PipelineError("FastQCJob did not complete successfully.")

    def _generate_job_script(self):
        lines = []

        details_file_name = f'{self.job_name}.array-details'
        sh_details_fp = join(self.output_path, details_file_name)

        lines.append("#!/bin/bash")

        job_name = f'{self.qiita_job_id}_{self.job_name}'
        lines.append(f"#SBATCH --job-name {job_name}")
        lines.append(f"#SBATCH -p {self.queue_name}")
        lines.append(f"#SBATCH -N {self.node_count}")
        lines.append(f"#SBATCH -n {self.nprocs}")
        lines.append("#SBATCH --time %d" % self.wall_time_limit)
        lines.append(f"#SBATCH --mem {self.jmem}")

        lines.append("#SBATCH --array 1-%d%%%d" % (
            len(self.commands), self.pool_size))

        lines.append("set -x")
        lines.append('date')
        lines.append('hostname')
        lines.append('echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}')
        lines.append(f'cd {self.output_path}')

        if self.modules_to_load:
            lines.append("module load " + ' '.join(self.modules_to_load))

        lines.append('offset=${SLURM_ARRAY_TASK_ID}')
        lines.append('step=$(( $offset - 0 ))')
        lines.append(f'cmd0=$(head -n $step {sh_details_fp} | tail -n 1)')
        lines.append('eval $cmd0')

        sentinel_file = f'{self.job_name}_$step.completed'
        lines.append(f'echo "Cmd Completed: $cmd0" > logs/{sentinel_file}')

        with open(self.job_script_path, 'w') as f:
            f.write('\n'.join(lines))

        with open(sh_details_fp, 'w') as f:
            f.write('\n'.join(self.commands))

    def parse_logs(self):
        log_path = join(self.output_path, 'logs')
        files = sorted(glob.glob(join(log_path, '*.out')))
        msgs = []

        for some_file in files:
            with open(some_file, 'r') as f:
                msgs += [line for line in f.readlines()
                         # note 'error' is not same
                         # requirement as found in QCJob.
                         # ('error:'). This is a very
                         # generalized filter.
                         if 'error' in line.lower()]

        return [msg.strip() for msg in msgs]
