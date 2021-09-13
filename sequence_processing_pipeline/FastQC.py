import os
from sequence_processing_pipeline.Job import Job
from os import makedirs, listdir
from os.path import join, dirname, abspath, exists, isfile, basename
from sequence_processing_pipeline.TorqueJob import TorqueJob
from time import sleep
from zipfile import ZipFile
import logging
import os
import re
import shutil


class FastQCJOb(Job):
    def __init__(self, root_dir, output_dir, nprocs, project):
        super().__init__()
        self.root_dir = root_dir
        self.output_dir = output_dir
        self.project = project
        self.nprocs = nprocs
        self.job_file_path = join(self.root_dir, 'Data/Fastq', self.project, 'fastqc_qsub.sh')

    def _generate_job_file(self, labname, chemistry, home, email_list, x, y):
        lines = []

        lines.append("#!/bin/bash")
        lines.append("#SBATCH --job-name=fastqc_%s_%s" % (self.project, os.path.basename(self.root_dir)))
        lines.append("#SBATCH --nodes=1")
        lines.append("#SBATCH --ntasks-per-node=%s" % self.nprocs)
        lines.append("#SBATCH --export=ALL")
        lines.append("#SBATCH --time=24:00:00")
        lines.append("#SBATCH --mem-per-cpu=6G")
        lines.append("#SBATCH --output=%s/fastqc_jobs/\%x-\%j.out" % (home, x, y))
        lines.append("#SBATCH --error=%s/fastqc_jobs/\%x-\%j.err" % (home, x, y))
        lines.append("#SBATCH --mail-type=ALL")
        lines.append("#SBATCH --mail-user=%s" % ','.join(email_list))
        # lines.append("source ~/miniconda3/bin/activate test_env_2")
        # lines.append("file=%s\%s" % (trim_file, slurm_array_task_id))
        lines.append('cmd="sh ~/seq_proc_dev/fastqc_parallel_bclconvert.sh %s %s %d %s %s %s"' % (
            self.root_dir, self.output_dir, self.nprocs, self.project, labname, chemistry))
        lines.append('echo "\$cmd"')
        lines.append("date")
        lines.append('echo "Executing: "\$cmd"')
        lines.append('eval "\$cmd"')
        lines.append("date")

        with open(self.job_file_path, 'w') as f:
            for line in lines:
                f.write("%s\n" % line)

    def run(self):
        # suffix=R*.fastq*
        # data_prefix=/projects/fastqc

        # fastq_raw=$fastq_output
        # atropos_qc_output=$fastq_output/*/atropos_qc
        # fastq_trimmed=$fastq_output/*/filtered_sequences

        self._generate_job_file(labname, chemistry, home, email_list, x, y)

        cmd = ['sbatch',
               '--parsable',
               '--qos=seq_proc',
               '--dependency=afterok:${pbs_job_id}',
               '--export=initial="true"',
               self.job_file_path]

        # TODO We may or may not want to wait
        fastqc_job_id = self.execute_sbatch_job_and_wait(cmd)

        # keep fastqc_parallel_bclconvert.sh as is for now.

