from os.path import join
from sequence_processing_pipeline.Job import Job
import logging
import os
import re


class FastQCJOb(Job):
    def __init__(self, root_dir, output_dir, nprocs, project, fastqc_path):
        super().__init__()
        self.root_dir = root_dir
        self.output_dir = output_dir
        self.job_script_path = join(self.root_dir, 'FASTQC.sh')
        self.project = project
        self.nprocs = nprocs
        self.stdout_log_path = join(self.root_dir, 'FASTQC_{}-{}.out.log')
        self.stderr_log_path = join(self.root_dir, 'FASTQC_{}-{}.err.log')
        self.fastqc_path = fastqc_path

    def _generate_job_script(self, queue_name, node_count, nprocs,
                             wall_time_limit, chemistry,  x, y):
        # filepath = ${seqdir}/Data/Fastq/${project}/fastqc_qsub.sh
        # TODO: Use Jinja template instead
        lines = []

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying
        # params w/environment variables is to do the same on QSUB but with
        # -v instead of --export. The syntax of the values are otherwise the
        # same. -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        job_name = "fastqc_%s_%s" % (
            self.project, os.path.basename(self.root_dir))
        lines.append("#PBS -N %s" % job_name)

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
        # using the larger value found in the two scripts (36 vs 24 hours)
        lines.append("#PBS -l walltime=%d:00:00" % wall_time_limit)

        # send an email to the list of users defined below when a job starts,
        # terminates, or aborts. This is used to confirm that the package's
        # own reporting mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # list of users to be contacted independently of this package's
        # notification system, when a job starts, terminates, or gets aborted.
        lines.append("#PBS -M ccowart@ucsd.edu,jdereus@ucsd.edu,qiita.help@gmail.com")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        # lines.append("#PBS -l pmem=10gb")

        # --output -> -o
        lines.append("#PBS -o %s" % self.stdout_log_path.format(x, y))
        lines.append("#PBS -e %s" % self.stderr_log_path.format(x, y))

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1

        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use root_dir instead.
        lines.append("set -x")
        lines.append("date '+%s' > /{}/{}.log".format(self.root_dir, job_name))
        lines.append("source ~/miniconda3/bin/activate test_env_2")
        lines.append("module load fastp_0.20.1 samtools_1.12 minimap2_2.18")
        lines.append("cd %s" % self.root_dir)
        # lines.append("export PATH=$PATH:/usr/local/bin")
        lines.append("#file=${trim_file}\${SLURM_ARRAY_TASK_ID}")
        lines.append('%s %s %s %s %s %s' % (
            self.fastqc_path, self.root_dir, self.output_dir, self.nprocs,
            self.project, chemistry))
        lines.append("echo $? >> /%s/%s.log" % (self.root_dir, job_name))
        lines.append("date '+%s' >> /{}/{}.log".format(
            self.root_dir, job_name))

        with open(self.job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % self.job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub('\s+', ' ', line)
                f.write("%s\n" % line)

    def run(self):
        pass


'''
#def run(self, queue_name, node_count, nprocs, wall_time_limit, chemistry, x,
y):
self._generate_job_script(queue_name, node_count, nprocs, wall_time_limit,
chemistry,  x, y)

job_info = self.qsub(self.job_script_path, None, None)
logging.info("Successful job: %s" % job_info)

parsable just means sbatch prints out the jobid and hte cluster name, kind of
like qsub does.
qos=sec_proc is just the quality of service profile. look for an  equivalent on
 our torque install.
depdendency=afterok:pbs_job_id means that this job is allowed to run only after
 pbs_job_id has completed successfully with a return code of 0.
    curiously, it doesn't say if the job queue aborts early if this dependency
     can never be met.
    it also says that you can't redo the job that failed successfully and have
    this job suddenly run. one time only
initial=true isn't needed, because initial is already reassigned inside the
script fastqc_qsub.sh calls.

fastqc_job_id =$(sbatch - -parsable - -qos
                 =seq_proc --dependency=afterok:${pbs_job_id}
                  --export=initial="true" ${seqdir} / Data / Fastq / ${project}
                  / fastqc_qsub.sh)

# TODO: if job returns successfully, notify user(s).
#  Users will be notified through PipelineErrors for
#  unsuccessful jobs.
'''
