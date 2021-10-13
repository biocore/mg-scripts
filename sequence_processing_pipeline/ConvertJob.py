from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os import walk
from os.path import join, abspath, exists, basename
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
import logging
import re


class ConvertJob(Job):
    def __init__(self, run_dir, sample_sheet_path, output_directory,
                 queue_name, node_count, nprocs, wall_time_limit, pmem,
                 bcl_tool_path, modules_to_load, qiita_job_id):
        '''
        ConvertJob provides a convenient way to run bcl-convert or bcl2fastq
        on a directory BCL files to generate Fastq files.
        :param run_dir: The 'run' directory that contains BCL files.
        :param sample_sheet_path: The path to a sample-sheet.
        :param output_directory: The path to store Fastq output.
        :param queue_name: The name of the Torque queue to use for processing.
        :param node_count: The number of nodes to request.
        :param nprocs: The maximum number of parallel processes to use.
        :param wall_time_limit: A hard time limit to bound processing.
        :param bcl_tool_path: The path to either bcl2fastq or bcl-convert.
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        '''
        self.job_name = 'ConvertJob'
        super().__init__(abspath(run_dir),
                         self.job_name,
                         [bcl_tool_path],
                         modules_to_load)
        self._directory_check(self.run_dir, create=False)
        self._validate_bcl_directory()
        self.modules_to_load = modules_to_load

        tmp = join(self.run_dir, 'Data', 'Fastq')
        # create the directory to store fastq files
        self._directory_check(tmp, create=True)
        # generate_job_script_path() also generates numbered logs for stdout
        # and stderr, but we won't need these. ConvertJob() just needs one
        # pair of logs.
        self.job_script_path, tmp1, tmp2 = self.generate_job_script_path()
        # self.run_id is of the form 210518_A00953_0305_AHCJT7DSX2
        self.run_id = basename(self.run_dir)
        self.job_name = f"{qiita_job_id}_ConvertJob"
        self.sample_sheet_path = sample_sheet_path
        self.output_directory = output_directory
        self.queue_name = queue_name
        self.node_count = node_count
        self.nprocs = nprocs
        self.wall_time_limit = wall_time_limit
        self.pmem = pmem
        self.bcl_tool = bcl_tool_path

        tmp = False
        for executable_name in ['bcl2fastq', 'bcl-convert']:
            if executable_name in self.bcl_tool:
                tmp = True
                break

        if not tmp:
            raise PipelineError(f'{self.bcl_tool} is not the path to a known'
                                'executable')

        self._file_check(self.sample_sheet_path)

        klss = KLSampleSheet(self.sample_sheet_path)

        if not validate_and_scrub_sample_sheet(klss):
            # if sample sheet is valid, we can assume it has the parameters
            # that we will need. No need to check for individual params.
            raise PipelineError(
                f"Sample sheet {self.sample_sheet_path} is not valid.")

        self._directory_check(self.output_directory, create=True)

        # required files for successful operation
        required_files = ['RTAComplete.txt', 'RunInfo.xml']
        for some_file in required_files:
            s = join(self.run_dir, some_file)
            self._file_check(s)

    def _validate_bcl_directory(self):
        '''

        :return:
        '''
        bcl_directory = join(self.run_dir, 'Data', 'Intensities', 'BaseCalls')
        bcl_count = 0
        if exists(bcl_directory):
            for root, dirs, files in walk(bcl_directory):
                for some_file in files:
                    # subdirectories and files other than '.bcl' files from
                    # bcl_directory are to be expected.
                    if some_file.endswith(('.bcl', '.cbcl')):
                        bcl_count += 1
        else:
            raise PipelineError(f"input_directory '{self.run_dir}' does not "
                                "contain subdirectory "
                                "Data/Intensities/BaseCalls.")

        if bcl_count == 0:
            raise PipelineError(f"input_directory '{bcl_directory}' does not "
                                f"contain bcl files ({bcl_count}).")

    def _generate_job_script(self):
        '''
        Generate a Torque job script for processing the supplied run_directory.
        :return: The path to the newly-created job-script.
        '''
        lines = []

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying
        # params w/environment variables is to do the same on QSUB but with -v
        # instead of --export. The syntax of the values are otherwise the same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        lines.append(f"#PBS -N {self.job_name}")

        # what torque calls a queue, slurm calls a partition
        # SBATCH -p SOMETHING -> PBS -q SOMETHING
        # (Torque doesn't appear to have a quality of service (-q) option. so
        # this will go unused in translation.)
        lines.append(f"#PBS -q {self.queue_name}")

        # request one node
        # Slurm --ntasks-per-node=<count> -> -l ppn=<count>	in Torque
        lines.append("#PBS -l nodes=%d:ppn=%d" % (self.node_count,
                                                  self.nprocs))

        # Slurm --export=ALL -> Torque's -V
        lines.append("#PBS -V")

        # Slurm
        # walltime limit --time=24:00:00 -> Torque's -l walltime=<hh:mm:ss>
        # using the larger value found in the two scripts (36 vs 24 hours)
        lines.append("#PBS -l walltime=%d:00:00" % self.wall_time_limit)

        # send an email to the list of users defined below when a job starts,
        # terminates, or aborts. This is used to confirm that the package's
        # own reporting mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # list of users to be contacted independently of this package's
        # notification system, when a job starts, terminates, or gets aborted.
        lines.append("#PBS -M qiita.help@gmail.com")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        lines.append(f"#PBS -l pmem={self.pmem}")

        # --output -> -o
        lines.append(f"#PBS -o {self.stdout_log_path}")
        lines.append(f"#PBS -e {self.stderr_log_path}")

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1
        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use run_dir instead.
        lines.append("set -x")
        lines.append(f'cd {self.run_dir}')
        tmp = "module load " + ' '.join(self.modules_to_load)
        lines.append(tmp)

        # Assume that the bcl-convert tool is named 'bcl-convert' and choose
        # accordingly.
        if 'bcl-convert' in self.bcl_tool:
            lines.append(('%s '
                          '--sample-sheet %s '
                          '--output-directory %s '
                          '--bcl-input-directory . '
                          '--bcl-num-decompression-threads 16 '
                          '--bcl-num-conversion-threads 16 '
                          '--bcl-num-compression-threads 16 '
                          '--bcl-num-parallel-tiles 16 '
                          '--bcl-sampleproject-subdirectories true '
                          '--force') % (self.bcl_tool,
                                        self.sample_sheet_path,
                                        self.output_directory))

            # equivalent cp for bcl-conversion (see below) needed.
        else:
            lines.append(('%s '
                          '--sample-sheet %s '
                          '--minimum-trimmed-read-length 1 '
                          '--mask-short-adapter-reads 1 '
                          '-R . '
                          '-o %s '
                          '--loading-threads 16 '
                          '--processing-threads 16 '
                          '--writing-threads 16 '
                          '--create-fastq-for-index-reads '
                          '--ignore-missing-positions ') %
                         (self.bcl_tool,
                          self.sample_sheet_path,
                          self.output_directory))

            # Stats directory is needed downstream.
            # It's a small directory (~64k)
            # Copy here instead of pulling it downstream.
            lines.append(f'cp -r ./Data/Fastq/Stats ./{self.run_id}')

        with open(self.job_script_path, 'w') as f:
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub(r'\s+', ' ', line)
                f.write(f"{line}\n")

        # for convenience
        return self.job_script_path

    def run(self):
        '''
        Run BCL2Fastq conversion on data in root directory.
        Job-related parameters are specified here for easy adjustment and
        job restart.
        :return:
        '''
        script_path = self._generate_job_script()
        logging.debug(script_path)

        job_info = self.qsub(self.job_script_path, None, None)
        logging.info(f'Successful job: {job_info}')
