from os import listdir
from os.path import join
from sequence_processing_pipeline.Job import Job
import logging


logging.basicConfig(level=logging.DEBUG)


class FastQCJob(Job):
    def __init__(self, run_dir, output_directory, nprocs, nthreads,
                 fastqc_path, modules_to_load, qiita_job_id, run_id,
                 queue_name, node_count, wall_time_limit, pmem):
        self.job_name = 'FastQCJob'
        super().__init__(run_dir, self.job_name, [fastqc_path],
                         modules_to_load)
        self.nprocs = nprocs
        self.nthreads = nthreads
        self.fastqc_path = fastqc_path
        self.modules_to_load = modules_to_load
        # assume that the path to raw fastqs is run_dir/Data/Fastq
        self.raw_fastq_path = join(self.run_dir, 'Data', 'Fastq')
        # output_directory is the root path for creating a FastQC directory
        # with FastQC output. It's also the location for finding trimmed and
        # filtered fastq files as well.
        self.fastqc_output_path = output_directory
        # '.../output_dir/run_id/project/filter_type'
        self.processed_fastq_path = join(self.run_dir, run_id)
        self.qiita_job_id = qiita_job_id
        self.trim_file = 'fqc_split_file_'

        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.products_dir = join(self.run_dir, run_id)
        self.pmem = pmem

    def get_cmds(self):
        results = []

        # gather the parameters for processing all relevant raw fastq files.
        params = self._scan_fastq_files(self.raw_fastq_path,
                                        self.fastqc_output_path,
                                        is_raw_input=True)

        for number_of_threads, file_path, output_path in params:
            cmd = ['fastqc', '--noextract', '-t', str(number_of_threads),
                   file_path, '-o', output_path]
            results.append(' '.join(cmd))

        # next, do the same for the trimmed/filtered fastq files.
        params = self._scan_fastq_files(self.processed_fastq_path,
                                        self.fastqc_output_path)

        for number_of_threads, file_path, output_path in params:
            cmd = ['fastqc', '--noextract', '-t', str(number_of_threads),
                   file_path, '-o', output_path]
            results.append(' '.join(cmd))

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
                    # results are a_trim = True, h_filter= = False
                    filter_type = 'trimmed_sequences'
                elif 'filtered_sequences' in tmp:
                    # results are a_trim = True, h_filter= = True
                    # (trimmed AND filtered)
                    filter_type = 'filtered_sequences'
                elif 'amplicon' in tmp:
                    # results are a_trim = False, h_filter= = False
                    filter_type = 'amplicon'
                else:
                    if is_raw_input:
                        filter_type = 'raw'
                    else:
                        raise ValueError("indeterminate type")

                if filter_type != 'amplicon':
                    # filter out index '_In_' files
                    files = [x for x in files if '_R1_' in x or '_R2_' in x]

                logging.debug("%s is a project directory" % project_dir)
                results.append((directory, filter_type, project_dir, files))

        return results

    def _scan_fastq_files(self, input_path, fastqc_output_path,
                          is_raw_input=False):
        '''
        Scan a .../run_id/Data/Fastq directory or an .../output_dir/run_id
        directory for FastQC job information.
        :param input_path: The path to a directory w/Fastq files.
        :param fastqc_output_path: The path to store all FastQC output.
        :param is_raw_input: Boolean to modify the output path.
        :param for_multiqc: Boolean. Generate params for FastQC or MultiQC.
        :return: A list of tuples suitable for adding to a subprocess.Pool
        '''
        projects = self._find_projects(input_path, is_raw_input)

        fastqc_results = []

        for project_name, filter_type, fastq_path, files in projects:
            tmp = join(fastqc_output_path, project_name)
            if is_raw_input:
                output_dir = join(tmp, 'bclconvert')
            else:
                output_dir = join(tmp, filter_type)
            # makedirs(output_dir, exist_ok=True)
            for some_file in files:
                fastqc_results.append((self.nthreads, some_file, output_dir))

        return fastqc_results

    def _generate_trim_files(self, fastq_files, split_count):
        def _chunk_list(some_list, n):
            # taken from https://bit.ly/38S9O8Z
            for i in range(0, len(some_list), n):
                yield some_list[i:i + n]

        # put these split files in the same location as the fastq_files for
        # the project. Assume all filepaths in fastq_files have the same
        # result for dirname().
        # destination_path = dirname(fastq_files[0])
        destination_path = self.run_dir

        new_files = []
        count = 0
        for chunk in _chunk_list(fastq_files, split_count):
            trim_file_name = '%s%d' % (self.trim_file, count)
            trim_file_path = join(destination_path, trim_file_name)
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
        cmds = self.get_cmds()
        split_count = self._generate_split_count(len(cmds))
        lines_per_split = int((len(cmds) + split_count - 1) / split_count)
        self._generate_trim_files(cmds, lines_per_split)

        script_path = self._generate_job_script(cmds,
                                                self.queue_name,
                                                self.node_count,
                                                self.nprocs,
                                                self.wall_time_limit,
                                                split_count,
                                                self.pmem)

        pbs_job_id = self.qsub(script_path, None, None)
        logging.debug(pbs_job_id)

    def _generate_job_script(self, cmds, queue_name, node_count, nprocs,
                             wall_time_limit, split_count, pmem):
        lines = []

        # unlike w/ConvertBCL2FastqJob, multiple qc.sh scripts will be
        # generated, one for each project defined in the sample sheet.
        # since we will generate a job-script for each project, we won't use
        # Job.stdout_log_path and stderr_log_path. Instead we'll use the
        # numbered ones provided by generate_job_script_path().
        job_script_path, output_log_path, error_log_path = \
            self.generate_job_script_path()

        # tf_fp = join(self.run_dir, self.trim_file + '0')
        sh_details_fp = join(self.run_dir, self.trim_file + '.array-details')

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying
        # params w/environment variables is to do the same on QSUB but with
        # -v instead of --export. The syntax of the values are otherwise the
        # same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        lines.append("#PBS -N %s" % f"{self.qiita_job_id}_FastQCJob")

        # what torque calls a queue, slurm calls a partition
        # SBATCH -p SOMETHING -> PBS -q SOMETHING
        # (Torque doesn't appear to have a quality of service (-q) option. so
        # this will go unused in translation.)
        lines.append("#PBS -q %s" % queue_name)

        # request one node
        # Slurm --ntasks-per-node=<count> -> -l ppn=<count>	in Torque
        lines.append("#PBS -l nodes=%d:ppn=%d" % (node_count, nprocs))

        # Slurm --export=ALL -> Torque's -V
        lines.append("#PBS -V")

        # Slurm
        # walltime limit --time=24:00:00 -> Torque's -l walltime=<hh:mm:ss>
        # using the larger value found in the two scripts (72 vs ? hours)
        lines.append("#PBS -l walltime=%d:00:00" % wall_time_limit)

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        lines.append(f"#PBS -l pmem={pmem}")

        # --output -> -o
        lines.append("#PBS -o localhost:%s.${PBS_ARRAYID}" % output_log_path)
        lines.append("#PBS -e localhost:%s.${PBS_ARRAYID}" % error_log_path)

        lines.append("#PBS -t 1-%d%%20" % len(cmds))

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1

        # We probably do not need to activate this Python environment, but I
        # will store it here in comments.
        # source ~/miniconda3/bin/activate test_env_2

        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use run_dir instead.
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
            f.write('\n'.join(cmds))

        return job_script_path
