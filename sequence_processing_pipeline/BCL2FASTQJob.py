import logging
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join, abspath, exists
import re
from os import makedirs, walk
from metapool import KLSampleSheet, validate_sample_sheet


class BCL2FASTQJob(Job):
    '''
    BCL2FASTQJob implements a way to run bcl2fastq on a directory of BCL
    files. It builds on TorqueJob's ability to push a job onto Torque and wait
    for it to finish.
    '''
    def __init__(self, root_dir, sample_sheet_path, output_directory, bcl2fastq_path):
        super().__init__()
        self.root_dir = abspath(root_dir)
        self._directory_check(self.root_dir, create=False)
        self._validate_bcl_directory()

        tmp = join(self.root_dir, 'Data', 'Fastq')
        makedirs(tmp, exist_ok=True)
        self.job_script_path = join(self.root_dir, 'BCL2FASTQ.sh')
        self.stdout_log_path = 'localhost:' + join(self.root_dir, 'BCL2FASTQ.out.log')
        self.stderr_log_path = 'localhost:' + join(self.root_dir, 'BCL2FASTQ.err.log')
        # self.unique_name is of the form 210518_A00953_0305_AHCJT7DSX2
        self.unique_name = self.root_dir.split('/')[-1]
        self.job_name = "BCL2FASTQ_%s" % self.unique_name
        self.sample_sheet_path = sample_sheet_path
        self.output_directory = output_directory
        self.bcl2fastq_path = bcl2fastq_path

        if not exists(self.sample_sheet_path):
            raise PipelineError("Sample sheet %s does not exist." % self.sample_sheet_path)

        klss = KLSampleSheet(self.sample_sheet_path)

        if not validate_sample_sheet(klss):
            # if sample sheet is valid, we can assume it has the parameters
            # that we will need. No need to check for individual params.
            raise PipelineError("Sample sheet %s is not valid." % self.sample_sheet_path)

        if not exists(self.output_directory):
            try:
                makedirs(self.output_directory, exist_ok=True)
            except OSError as e:
                raise PipelineError("Cannot create output directory %s. %s" % (self.output_directory, str(e)))

        if not exists(self.bcl2fastq_path):
            raise PipelineError("Path to bcl2fastq binary '%s' is not valid." % self.bcl2fastq_path)

    def _validate_bcl_directory(self):
        bcl_directory = join(self.root_dir, 'Data', 'Intensities', 'BaseCalls')
        bcl_count = 0
        if exists(bcl_directory):
            for root, dirs, files in walk(bcl_directory):
                for some_file in files:
                    # subdirectories and files other than '.bcl' files from
                    # bcl_directory are to be expected.
                    if some_file.endswith('.bcl') or some_file.endswith('.cbcl'):
                        bcl_count += 1
        else:
            logging.debug("input_directory '%s' does not contain subdirectory Data/Intensities/BaseCalls." % self.root_dir)
            raise PipelineError("input_directory '%s' does not contain subdirectory Data/Intensities/BaseCalls." % self.root_dir)

        if bcl_count < 1:
            # we can increase to whatever threshold is acceptable.
            raise PipelineError("input_directory '%s' does not contain enough bcl files (%d)." % (bcl_directory, bcl_count))

    def _generate_job_script(self, queue_name, node_count, nprocs, wall_time_limit):
        # TODO: Use Jinja template instead
        lines = []

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying params
        # w/environment variables is to do the same on QSUB but with -v instead of
        # --export. The syntax of the values are otherwise the same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        lines.append("#PBS -N %s" % self.job_name)

        # what torque calls a queue, slurm calls a partition
        # SBATCH -p SOMETHING -> PBS -q SOMETHING
        # (Torque doesn't appear to have a quality of service (-q) option. so this
        # will go unused in translation.)
        lines.append("#PBS -q %s" % queue_name)

        # request one node
        # Slurm --ntasks-per-node=<count> -> -l ppn=<count>	in Torque
        lines.append("#PBS -l nodes=%d:ppn=%d" % (node_count, nprocs))

        # Slurm --export=ALL -> Torque's -V
        lines.append("#PBS -V")

        # Slurm walltime limit --time=24:00:00 -> Torque's -l walltime=<hh:mm:ss>
        # using the larger value found in the two scripts (36 vs 24 hours)
        lines.append("#PBS -l walltime=%d:00:00" % wall_time_limit)

        # send email to charlie when a job starts and when it terminates or
        # aborts. This is used to confirm the package's own reporting
        # mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # specify your email address
        # TODO: Send an email to jeff, and qiita.help as well.
        lines.append("#PBS -M ccowart@ucsd.edu")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        #lines.append("#PBS -l pmem=10gb")

        # --output -> -o
        lines.append("#PBS -o %s" % self.stdout_log_path)
        lines.append("#PBS -e %s" % self.stderr_log_path)

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1
        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use root_dir instead.
        lines.append("set -x")
        lines.append("cd %s" % self.root_dir)
        lines.append("date '+%s' > bcl2fastq.job.log")
        lines.append("module load bcl2fastq_2.20.0.422")
        lines.append('%s --sample-sheet %s \
                      --mask-short-adapter-reads 1 \
                      -R . \
                      -o %s \
                      --loading-threads 8 \
                      --processing-threads 8 \
                      --minimum-trimmed-read-length 1 \
                      --writing-threads 2 \
                      --create-fastq-for-index-reads \
                      --ignore-missing-bcls' % (self.bcl2fastq_path,
                                                 self.sample_sheet_path,
                                                 self.output_directory))
        lines.append("echo $? >> bcl2fastq.job.log")
        lines.append("date '+%s' >> bcl2fastq.job.log")

        with open(self.job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % self.job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub('\s+', ' ', line)
                f.write("%s\n" % line)

    def run(self, queue_name, node_count, nprocs, wall_time_limit):
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
        self._generate_job_script(queue_name, node_count, nprocs, wall_time_limit)

        job_info = self.qsub(self.job_script_path, None, None)
        logging.info("Successful job: %s" % job_info)

        # TODO: if job returns successfully, notify user(s).
        #  Users will be notified through PipelineErrors for
        #  unsuccessful jobs.
