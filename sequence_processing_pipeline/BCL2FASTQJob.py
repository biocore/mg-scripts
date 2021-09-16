import logging
from sequence_processing_pipeline.TorqueJob import TorqueJob
from os.path import join, abspath
import re
from os import makedirs


class BCL2FASTQJob(TorqueJob):
    '''
    BCL2FASTQJob implements a way to run bcl2fastq on a directory of BCL
    files. It builds on TorqueJob's ability to push a job onto Torque and wait
    for it to finish.
    '''
    def __init__(self, root_dir, sample_sheet_path, output_directory, bcl2fastq_path):
        super().__init__()
        self.root_dir = abspath(root_dir)
        tmp = join(self.root_dir, 'Data', 'Fastq')
        makedirs(tmp, exist_ok=True)
        self.job_script_path = join(self.root_dir, 'BCL2FASTQ.sh')
        self.stdout_log_path = join(self.root_dir, 'BCL2FASTQ.out.log')
        self.stderr_log_path = join(self.root_dir, 'BCL2FASTQ.err.log')
        # self.unique_name is of the form 210518_A00953_0305_AHCJT7DSX2
        self.unique_name = self.root_dir.split('/')[-1]
        self.job_name = "BCL2FASTQ_%s" % self.unique_name
        self.sample_sheet_path = sample_sheet_path
        self.output_directory = output_directory
        self.bcl2fastq_path = bcl2fastq_path

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
        lines.append("date '+%s' > /{}/{}.log".format(self.root_dir, self.job_name))
        lines.append("module load bcl2fastq_2.20.0.422")
        lines.append("cd %s"  % self.root_dir)
        # lines.append("export PATH=$PATH:/usr/local/bin")
        lines.append('%s --sample-sheet %s --mask-short-adapter-reads \
                      1 -R . -o %s --loading-threads 8 --processing-threads \
                      --minimum-trimmed-read-length 1 \
                      8 --writing-threads 2 --create-fastq-for-index-reads \
                      --ignore-missing-bcls' % (self.bcl2fastq_path,
                                                 self.sample_sheet_path,
                                                 self.output_directory))
        lines.append("echo $? >> /%s/%s.log" % (self.root_dir, self.job_name))
        lines.append("date '+%s' >> /{}/{}.log".format(self.root_dir, self.job_name))

        with open(self.job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % self.job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub('\s+', ' ', line)
                f.write("%s\n" % line)

    def run(self, queue_name, node_count, nprocs, wall_time_limit):
        self._generate_job_script(queue_name, node_count, nprocs, wall_time_limit)

        job_info = self.qsub(self.job_script_path, None, None)
        logging.info("Successful job: %s" % job_info)

        # TODO: if job returns successfully, notify user(s).
        #  Users will be notified through PipelineErrors for
        #  unsuccessful jobs.
