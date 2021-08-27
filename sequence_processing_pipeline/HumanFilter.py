import os

class HumanFilter:
    def __init__(self, seq_dir):
        self.seq_dir = seq_dir

    def qsub_file_generator(self, seq_dir, project, file_base, nprocs, home, email_list, trim_file, slurm_array_task_id,
                            chemistry, output_dir, base_name, adapter_a, adapter_A, a_trim, h_filter, qiita_proj,
                            final_output_dir, cmd):
        lines = []

        # TODO: replace with a stock file
        lines.append("#!/bin/bash")
        lines.append("#SBATCH --job-name=%s_\%A_\%a" % project)
        lines.append("#SBATCH --ntasks=%d" % nprocs)
        lines.append("#SBATCH --ntasks-per-node=%d" % nprocs)
        lines.append("#SBATCH --export=ALL")
        # TODO: assuming for now that 72 hours for all jobs is appropriate
        lines.append("#SBATCH --time=72:00:00")
        # TODO: ditto with other parameters
        lines.append("#SBATCH --mem-per-cpu=6G")
        lines.append("#SBATCH --output=%s/filter_jobs/\%x_\%A_\%a.out" % home)
        lines.append("#SBATCH --error=%s/filter_jobs/\%x_\%A_\%a.err" % home)
        lines.append("#SBATCH --mail-type=ALL")
        lines.append("#SBATCH --mail-user=%s" % ','.join(email_list))
        # TODO: assume for now this hardcoded path is acceptable
        lines.append("source ~/miniconda3/bin/activate test_env_2")
        # TODO: assume for now the backslash below is intentional
        lines.append("file=%s\%s" % (trim_file, slurm_array_task_id))
        lines.append(
            "cmd=\"sh ~/seq_proc_dev/fpmmp.sh -d %s -D /Data/Fastq -S %s\%s -p %s -C %s -c %d -o %s -O %s -a %s -A %s -g %s -G %s -q %s -f %s\"" % (
            seq_dir, trim_file, slurm_array_task_id, project, chemistry, nprocs, output_dir, base_name,
            os.path.basename(seq_dir), adapter_a, adapter_A, a_trim, h_filter, qiita_proj, final_output_dir))
        lines.append("echo \"%s\"" % cmd)
        lines.append("date")
        lines.append("echo 'Executing: %s'" % cmd)
        lines.append("eval \"%s\"" % cmd)
        lines.append("date")

        file_path = os.path.join(seq_dir, 'Data/Fastq', project, file_base + '_qsub.sh')
        with open(file_path, 'w') as f:
            for line in lines:
                f.write("%s\n" % line)

    def human_filter(self, seq_dir, output_dir, fastq_output):
        '''
        Called with parameters: human_filter "${seqdir}" "$output_dir" "${fastq_output}"
        :param self:
        :param seq_dir:
        :param output_dir:
        :param fastq_output:
        :return:
        '''
        # if seq_dir doesn't contain a run_config.txt file, then abort
        # it looks like fastq_output might be overwritten to be seq_dir/Data/Fastq

        # check run_config.txt for duplicate line_per_split
        # possible same project run on multiple samplesheets
        run_config_path = os.path.join(seq_dir, 'run_config.txt')
        lines = None
        with open(run_config_path, 'r') as f:
            lines = f.readlines()
            lines = list(set([x.strip() for x in lines]))
        if lines:
            with open(run_config_path, 'w') as f:
                for line in lines:
                    f.write('%s\n' % line)
        else:
            raise ValueError("run_config.txt file could not be processed.")

        # SECTIONS COMMENTED OUT NEED TO BE PORTED, NOT REMOVED!

        '''
        while IFS=" " read project adapter_a adapter_A a_trim h_filter qiita_proj; do

        proj_check=$(echo $project | awk -F_ '{print $NF}')
        echo $project $adapter_a $adapter_A $a_trim $h_filter $qiita_proj

        ### remove any trailing newlines (^M) from variable
        qiita_proj=$(echo -n $qiita_proj)
        qiita_proj=$(echo "$qiita_proj"|tr -d '\r')
        pushd ${fastq_output}/${project}
        input_count=$(ls ${fastq_output}/${project}/*R1*.fastq.gz | wc -l)
        '''
        input_count = 1
        trim_file = 1

        if input_count > 2000:
            split_count = 16
        elif input_count <= 2000 and input_count > 1000:
            split_count = 10
        elif input_count <= 1000 and input_count > 500:
            split_count = 4
        else:
            split_count = 1

        if trim_file == None or trim_file == '':
            trim_file = 'split_file_'  # or possibly the variable split_file_

        job_count = split_count + 1

        file_base = os.path.basename(seq_dir)

        for root, dirs, files in os.walk(whatever_directory_trim_files_are_in):
            for some_file in files:
                some_path = os.path.join(root, some_file)
                if some_path.startswith(trim_file):
                    os.remove(some_path)
                    # TODO: Note Jeff's rm -fv (all trim files) will be faster than this

        '''
        find . -maxdepth 1 -name "*.fastq.gz" -type f | grep "_R1_" | cut -f2 -d"/" > ${file_base}_file_list.txt
          #ls | grep *.fastq.gz > ${file_base}_file_list.txt
          line_count=$(cat ${file_base}_file_list.txt | wc -l)
          line_per_split=$(( $(( $line_count + $split_count -1 ))/$split_count ))

          split -l $line_per_split -d ${file_base}_file_list.txt split_file_; rename ${trim_file}0 ${trim_file} ${trim_file}*
          if [[ ${line_count} -eq "0" ]]; then
            :
          fi
        '''

        self.qsub_file_generator(seq_dir, project, file_base, nprocs, home, email_list, trim_file, slurm_array_task_id,
                                 chemistry, output_dir, base_name, adapter_a, adapter_A, a_trim, h_filter, qiita_proj,
                                 final_output_dir, cmd)

        '''

        pbs_job_id=$(sbatch --qos=seq_proc --parsable --array=0-$(($split_count - 1)) ${seqdir}/Data/Fastq/${project}/${file_base}_qsub.sh)

        ### this should never fail as it only is executed while there is another project to process
        ### first time through
        #if [[ -z $pbs_job_list ]]; then
          #pbs_job_list="-W depend=afterokarray:$pbs_job_id"
          pbs_job_array="--dependency=afterok:$pbs_job_id"

        if [[ -z proj_array ]]; then
          proj_array=$project
        else
          proj_array+=$project
        fi

        ### submit fastqc for current project directory
        fastqc_process "${seqdir}" "${output_dir}" "${fastq_output}" "${pbs_job_id}" "${pbs_job_array}" "${proj_array[@]}"

        ### out of project directory
        popd

        done < ${seqdir}/uniq_run_config.txt
    	popd
        '''







