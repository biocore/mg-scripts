


def fastqc_process_file_generator(self, seq_dir, project, file_base, nprocs, home, email_list, trim_file, slurm_array_task_id,
                        chemistry, output_dir, base_name, adapter_a, adapter_A, a_trim, h_filter, qiita_proj,
                        final_output_dir, cmd):
    # output filepath is ${seqdir}/Data/Fastq/${project}/fastqc_qsub.sh

    lines = []

    # TODO: replace with a stock file
    lines.append("#!/bin/bash")
    lines.append("#SBATCH --job-name=fastqc_%s_$(basename ${%s})" % (project, seq_dir))
    lines.append("#SBATCH --nodes=1")
    lines.append("#SBATCH --ntasks-per-node=%s" % nprocs)
    lines.append("#SBATCH --export=ALL")
    lines.append("#SBATCH --time=24:00:00")
    lines.append("#SBATCH --mem-per-cpu=6G")
    lines.append("#SBATCH --output=%s/fastqc_jobs/%x-%j.out" % (home))
    lines.append("#SBATCH --error=%s/fastqc_jobs/%x-%j.err" % (home))
    lines.append("#SBATCH --mail-type=ALL")
    lines.append("#SBATCH --mail-user=%s" % ','.join(email_list))
    lines.append("source ~/miniconda3/bin/activate test_env_2")
    lines.append("file=%s\%s" % (trim_file, slurm_array_task_id))
    lines.append(
        "cmd=\"sh ~/seq_proc_dev/fastqc_parallel_bclconvert.sh ${seqdir} ${output_dir} ${NPROCS} ${project} ${labname} ${chemistry}\"")
    lines.append("echo \"\$cmd\"")
    lines.append("date")
    lines.append("echo \"Executing: \"$cmd\" \"")
    lines.append("eval \"\$cmd\"")
    lines.append("date")
    lines.append("#( printf '%s\n\n' \"fastqc processing complete for ${project}_$(basename ${seqdir})\"")
    lines.append(
        "#  printf '%s\n' \"fastqc results located at http://kl-fastqc.ucsd.edu/fastqc/${labname}/$(basename ${seqdir})\"")
    lines.append("#) | mailx -s fastqc complete -r jdereus@barnacle.ucsd.edu $fastqc_mail")

    file_path = os.path.join(seq_dir, 'Data/Fastq', project, file_base + '_qsub.sh')
    with open(file_path, 'w') as f:
        for line in lines:
            f.write("%s\n" % line)


'''
function fastqc_process() {
  suffix=R*.fastq*
  data_prefix=/projects/fastqc

  fastq_raw=$fastq_output
  atropos_qc_output=$fastq_output/*/atropos_qc
  fastq_trimmed=$fastq_output/*/filtered_sequences

###for proj in ${proj_array[@]}; do
  cat > ${seqdir}/Data/Fastq/${project}/fastqc_qsub.sh <<-EOF
#!/bin/bash
#SBATCH --job-name=fastqc_${project}_$(basename ${seqdir})
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${NPROCS}
#SBATCH --export=ALL
#SBATCH --time=24:00:00
###96:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --output=${HOME}/fastqc_jobs/%x-%j.out
#SBATCH --error=${HOME}/fastqc_jobs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${email_final[@]}
###jdereus@ucsd.edu
####,${maillist}

source ~/miniconda3/bin/activate test_env_2

file=${trim_file}\${SLURM_ARRAY_TASK_ID}

NOTE: fastqc_parallel_bclconvert.sh needed

cmd="sh ~/seq_proc_dev/fastqc_parallel_bclconvert.sh ${seqdir} ${output_dir} ${NPROCS} ${project} ${labname} ${chemistry}"
echo "\$cmd"
date
echo "Executing: \"$cmd\" "
eval "\$cmd"
date

#( printf '%s\n\n' "fastqc processing complete for ${project}_$(basename ${seqdir})"
#  printf '%s\n' "fastqc results located at http://kl-fastqc.ucsd.edu/fastqc/${labname}/$(basename ${seqdir})"
#) | mailx -s "fastqc complete" -r jdereus@barnacle.ucsd.edu $fastqc_mail
EOF

  fastqc_job_id=$(sbatch --parsable --qos=seq_proc --dependency=afterok:${pbs_job_id} --export=initial="true" ${seqdir}/Data/Fastq/${project}/fastqc_qsub.sh)
}
'''