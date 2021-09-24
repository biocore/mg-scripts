import logging
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join, abspath, exists
import re
from os import makedirs, walk
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet


class ConvBCL2FastqJob(Job):
    def __init__(self, run_dir, sample_sheet_path, output_directory,
                 bcl_executable_path, use_bcl_convert, queue_name, node_count,
                 nprocs, wall_time_limit, complete_runs_path=None):
        '''
        ConvBCL2FastqJob implements a way to run bcl2fastq on a directory of
        BCL files.
        Submit a Torque job where the contents of run_dir are processed using
        fastp, minimap, and samtools. Human-genome sequences will be filtered
        out if needed.
        :param run_dir:
        :param sample_sheet_path:
        :param output_directory:
        :param bcl_executable_path:
        :param use_bcl_convert:
        :param queue_name:
        :param node_count:
        :param nprocs:
        :param wall_time_limit:
        :param complete_runs_path:
        '''
        super().__init__()
        self.run_dir = abspath(run_dir)
        self._directory_check(self.run_dir, create=False)
        self._validate_bcl_directory()

        tmp = join(self.run_dir, 'Data', 'Fastq')
        makedirs(tmp, exist_ok=True)
        self.job_script_path = join(self.run_dir, 'ConvertBCL2Fastq.sh')
        self.stdout_log_path = 'localhost:' + join(
            self.run_dir, 'ConvertBCL2Fastq.out.log')
        self.stderr_log_path = 'localhost:' + join(
            self.run_dir, 'ConvertBCL2Fastq.err.log')
        # self.unique_name is of the form 210518_A00953_0305_AHCJT7DSX2
        self.unique_name = self.run_dir.split('/')[-1]
        self.job_name = "ConvertBCL2Fastq_%s" % self.unique_name
        self.sample_sheet_path = sample_sheet_path
        self.output_directory = output_directory
        self.bcl_executable_path = bcl_executable_path
        self.use_bcl_convert = use_bcl_convert
        self.queue_name = queue_name
        self.node_count = node_count
        self.nprocs = nprocs
        self.wall_time_limit = wall_time_limit
        self.complete_runs_path = complete_runs_path

        d = {'bcl2fastq': {}, 'bcl-convert': {}}
        d['bcl2fastq']['module_load'] = 'module load bcl2fastq_2.20.0.422'
        d['bcl-convert']['module_load'] = 'module load bclconvert_3.7.5'
        d['bcl2fastq']['executable_path'] = self.bcl_executable_path
        d['bcl-convert']['executable_path'] = self.bcl_executable_path
        d['bcl2fastq']['queue_name'] = ''
        d['bcl-convert']['queue_name'] = 'highmem'
        d['bcl2fastq']['output_log_name'] = ('localhost: %s/BCL2FASTQ.out.log'
                                             % self.run_dir)
        d['bcl2fastq']['error_log_name'] = ('localhost: %s/BCL2FASTQ.err.log'
                                            % self.run_dir)
        tmp = 'localhost: %s/BCLCONVERT.out.log' % self.run_dir
        d['bcl-convert']['output_log_name'] = tmp
        tmp = 'localhost: %s/BCLCONVERT.err.log' % self.run_dir
        d['bcl-convert']['error_log_name'] = tmp

        if use_bcl_convert:
            self.bcl_tool = d['bcl-convert']
        else:
            self.bcl_tool = d['bcl2fastq']

        self._file_check(self.sample_sheet_path)

        klss = KLSampleSheet(self.sample_sheet_path)

        if not validate_and_scrub_sample_sheet(klss):
            # if sample sheet is valid, we can assume it has the parameters
            # that we will need. No need to check for individual params.
            raise PipelineError(
                "Sample sheet %s is not valid." % self.sample_sheet_path)

        self._directory_check(self.output_directory, create=True)

        self._file_check(self.bcl_executable_path)

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
            raise PipelineError("input_directory '%s' does not contain "
                                "bcl files (%d)." % (bcl_directory, bcl_count))

    def _generate_job_script(self):
        '''

        :return:
        '''
        lines = []

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying
        # params w/environment variables is to do the same on QSUB but with -v
        # instead of --export. The syntax of the values are otherwise the same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        lines.append("#PBS -N %s" % self.job_name)

        # what torque calls a queue, slurm calls a partition
        # SBATCH -p SOMETHING -> PBS -q SOMETHING
        # (Torque doesn't appear to have a quality of service (-q) option. so
        # this will go unused in translation.)
        lines.append("#PBS -q %s" % self.queue_name)

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
        lines.append("#PBS -M "
                     "ccowart@ucsd.edu,"
                     "jdereus@ucsd.edu,"
                     "qiita.help@gmail.com")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        lines.append("#PBS -l pmem=10gb")

        # --output -> -o
        lines.append("#PBS -o %s" % self.bcl_tool['output_log_name'])
        lines.append("#PBS -e %s" % self.bcl_tool['error_log_name'])

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1
        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use run_dir instead.
        lines.append("set -x")
        lines.append("cd %s" % self.run_dir)
        lines.append(self.bcl_tool['module_load'])

        if self.use_bcl_convert:
            lines.append('%s --sample-sheet %s \
                          --mask-short-adapter-reads 1 \
                          -R . \
                          -o %s \
                          --loading-threads 8 \
                          --processing-threads 8 \
                          --minimum-trimmed-read-length 1 \
                          --writing-threads 2 \
                          --create-fastq-for-index-reads \
                          --ignore-missing-bcls' % (
                self.bcl_tool['executable_path'], self.sample_sheet_path,
                self.output_directory))
        else:
            lines.append('%s \
                          --sample-sheet %s \
                          --output-directory %s \
                          --bcl-input-directory . \
                          --bcl-num-decompression-threads 8 \
                          --bcl-num-conversion-threads 8 \
                          --bcl-num-compression-threads 16 \
                          --bcl-num-parallel-tiles 8 \
                          --bcl-sampleproject-subdirectories true --force' % (
                self.bcl_tool['executable_path'], self.sample_sheet_path,
                self.output_directory))

        if self.complete_runs_path:
            stats_path = join(self.output_directory, 'Stats', 'Stats.json')
            lines.append('cp %s %s' % (stats_path, self.complete_runs_path))

        with open(self.job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % self.job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub(r'\s+', ' ', line)
                f.write("%s\n" % line)

    def run(self):
        '''
        Run BCL2Fastq conversion on data in root directory.
        Job-related parameters are specified here for easy adjustment and
        job restart.
        :param queue_name:
        :param node_count:
        :param nprocs:
        :param wall_time_limit:
        :return:
        '''
        self._generate_job_script()

        job_info = self.qsub(self.job_script_path, None, None)
        logging.info("Successful job: %s" % job_info)

        # TODO: if job returns successfully, notify user(s).
        #  Users will be notified through PipelineErrors for
        #  unsuccessful jobs.
