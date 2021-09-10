import logging
from sequence_processing_pipeline.TorqueJob import TorqueJob
from os.path import join, basename, dirname, abspath
from os import makedirs
import re


class BCL2FASTQJob(TorqueJob):
    '''
    BCL2FASTQJob implements a way run bcl2fastq on a directory of BCL files.
    It builds on TorqueJob's ability to push a job onto Torque and wait for it
     to finish.
    '''
    def __init__(self, root_dir, project, sample_sheet_path, output_directory):
        super().__init__()
        logging.debug("BCL2FASTQJob Constructor called")
        self.root_dir = abspath(root_dir)
        self.project = project
        self.job_script_path = join(self.root_dir, 'Data', 'Fastq', project, 'fastqc_qsub.sh')
        self.job_name = "fastqc_%s_%s" % (project, basename(self.root_dir))
        self.nprocs = 16
        self.sample_sheet_path = sample_sheet_path
        self.stdout_log_path = join(dirname(self.job_script_path), 'BCL2FASTQ.out.log')
        self.stderr_log_path = join(dirname(self.job_script_path), 'BCL2FASTQ.err.log')
        self.output_dir_path = abspath(output_directory)
        self.bcl2fastq_bin_path = "~/bcl2fastq_2_20/bin/bcl2fastq"
        makedirs(dirname(self.job_script_path), exist_ok=True)
        makedirs(dirname(self.output_dir_path), exist_ok=True)

    def _make_job_script(self):
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
        lines.append("#PBS -q long8gb")

        # request one node
        # Slurm --ntasks-per-node=<count> -> -l ppn=<count>	in Torque
        lines.append("#PBS -l nodes=1:ppn=%d" % self.nprocs)

        # Slurm --export=ALL -> Torque's -V
        lines.append("#PBS -V")

        # Slurm walltime limit --time=24:00:00 -> Torque's -l walltime=<hh:mm:ss>
        # using the larger value found in the two scripts (36 vs 24 hours)
        lines.append("#PBS -l walltime=36:00:00")

        # send email to charlie when a job starts and when it terminates or
        # aborts. This is used to confirm the package's own reporting
        # mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # specify your email address
        lines.append("#PBS -M ccowart@ucsd.edu")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        #lines.append("#PBS -l pmem=10gb")

        # --output -> -o
        lines.append("#PBS -o %s" % self.stdout_log_path)
        lines.append("#PBS -e %s" % self.stderr_log_path)

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1

        # We probably do not need to activate this Python environment, but I
        # will store it here in comments.
        # source ~/miniconda3/bin/activate test_env_2

        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use root_dir instead.
        lines.append("cd %s"  % self.root_dir)
        lines.append('cmd="%s --sample-sheet %s --mask-short-adapter-reads \
                      1 -R . -o %s --loading-threads 8 --processing-threads \
                      8 --writing-threads 2 --create-fastq-for-index-reads \
                      --ignore-missing-bcls"' % (self.bcl2fastq_bin_path,
                                                 self.sample_sheet_path,
                                                 self.output_dir_path))
        lines.append("echo $cmd")
        lines.append("eval $cmd")
        lines.append("return_code=$?")
        p = "%s/%s.return_code" % (self.root_dir, self.job_name)
        lines.append("echo $returncode > %s" % p)
        lines.append("date >> %s" % p)

        with open(self.job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % self.job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub('\s+', ' ', line)
                f.write("%s\n" % line)

    def run(self):
        # we may actually need to make the job script in a sub-dir.
        self._make_job_script()
        job_info = self.qsub(self.job_script_path, None, None)

        # if the job returned successfully, we may want to send an
        # email to the user. If an error occurs, the user will
        # be notified by the email triggered by PipelineError().


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    job = BCL2FASTQJob('./good-bcl-directory',
                       'THDMI_US_99999',
                       './good-bcl-directory/good-sample-sheet.csv',
                       './good-bcl-directory/Data/Fastq/Output')
    job.run()


