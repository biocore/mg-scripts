import os
from sequence_processing_pipeline.Job import Job


class FastQCJOb(Job):
    def __init__(self):
        super().__init__()

    def _create_seq_proc_file(self, seq_dir, project, file_base, nprocs, labname, chemistry, home, email_list,
                              trim_file, slurm_array_task_id):
        file_path = os.path.join(seq_dir, 'Data', 'Fastq', project, 'fastqc_qsub.sh')

        lines = []

        with open(file_path, 'w') as f:
            lines.append("#!/bin/bash")
            lines.append("#SBATCH --job-name=fastqc_%s_%s" % (project, os.path.basename(seq_dir)))
            lines.append("#SBATCH --nodes=1")
            lines.append("#SBATCH --ntasks-per-node=%s" % nprocs)
            lines.append("#SBATCH --export=ALL")
            lines.append("#SBATCH --time=24:00:00")
            lines.append("#SBATCH --mem-per-cpu=6G")
            lines.append("#SBATCH --output=%s/fastqc_jobs/\%x-\%j.out" % (home))
            lines.append("#SBATCH --error=%s/fastqc_jobs/\%x-\%j.err" % (home))
            lines.append("#SBATCH --mail-type=ALL")
            lines.append("#SBATCH --mail-user=%s" % ','.join(email_list))
            lines.append("source ~/miniconda3/bin/activate test_env_2")
            lines.append("file=%s\%s" % (trim_file, slurm_array_task_id))
            lines.append('cmd="sh ~/seq_proc_dev/fastqc_parallel_bclconvert.sh %s %s %d %s %s %s"' % (
            seq_dir, output_dir, nprocs, project, labname, chemistry))
            lines.append('echo "\$cmd"')
            lines.append("date")
            lines.append('echo "Executing: "\$cmd"')
            lines.append('eval "\$cmd"')
            lines.append("date")

            # these were commented out and I don't think they're needed at this time.
            lines.append("#( printf '%s\n\n' \"fastqc processing complete for ${project}_$(basename ${seqdir})\"")
            lines.append(
                "#  printf '%s\n' \"fastqc results located at http://kl-fastqc.ucsd.edu/fastqc/${labname}/$(basename ${seqdir})\"")
            lines.append("#) | mailx -s fastqc complete -r jdereus@barnacle.ucsd.edu $fastqc_mail")

            file_path2 = os.path.join(seq_dir, 'Data/Fastq', project, 'fastqc_qsub.sh')

            with open(file_path2, 'w') as f:
                for line in lines:
                    f.write("%s\n" % line)

    def run(self):
        # suffix=R*.fastq*
        # data_prefix=/projects/fastqc

        # fastq_raw=$fastq_output
        # atropos_qc_output=$fastq_output/*/atropos_qc
        # fastq_trimmed=$fastq_output/*/filtered_sequences

        self._create_seq_proc_file()

        cmd = ['sbatch',
               '--parsable',
               '--qos=seq_proc',
               '--dependency=afterok:${pbs_job_id}',
               '--export=initial="true"',
               '${seqdir}/Data/Fastq/${project}/fastqc_qsub.sh']

        # TODO We may or may not want to wait
        fastqc_job_id = self.execute_sbatch_job_and_wait(cmd)

        # keep fastqc_parallel_bclconvert.sh as is for now.
