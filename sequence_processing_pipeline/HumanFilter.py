import os
import logging



class HumanFilter:
    def __init__(self, SequenceDirectoryObject):
        self.sdo = SequenceDirectoryObject
        self.seq_dir = self.sdo.seq_dir

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

        ### submit fastqc for current project directory
        # i think the issue here is that human filtering gets sbatched and then fastqc is run on the results?
        # so if sbatch is just queueing the job that's fine but it seems like fastqc relies on the output of
        # that so they're tied together.
        fastqc_process_file_generator "${seqdir}" "${output_dir}" "${fastq_output}" "${pbs_job_id}" "${pbs_job_array}" "${proj_array[@]}"
    	


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

class HumanFilter2:
    def __init__(self, seq_dir):
        self.csv_files = []
        self.seq_dir = seq_dir
        self.run_config_file_name = 'run_config.txt'
        self.run_config_file_path = os.path.join(self.seq_dir, self.run_config_file_name)

    def _run_param(self):
        skip_these_files = ['sav.csv']

        for csv_file in self.csv_files:
            if csv_file in skip_these_files:
                logging.debug("Skipping %s..." % csv_file)
            else:
                with open(csv_file, 'r') as f:
                    # we can (and should) open up this file as a proper INI
                    # file. However it may have trailing ',' characters,
                    # which may impede this. Hence, open it as a regular
                    # file and clean it up before processing.
                    lines = f.readlines()
                    # first, let's strip out \n and \r\n and any leading or
                    # trailing whitespaces
                    lines = [x.strip() for x in lines]
                    # second, let's strip out trailing ',' characters, if
                    # they're present. We'll keep any leading ',' characters
                    # or in-line ',' characters.
                    lines = [x.rstrip(',') for x in lines]
                    # lastly, there have been cases where lines contain only
                    # ',' characters. Hence, look for and remove these now
                    # empty lines before continuing.
                    lines = [x for x in lines if x]
                    # since the file is already in memory, we won't use INI
                    # library to parse it, (for now). Although it would be
                    # cleaner.
                    metadata = []
                    sentinel = False
                    for i in range(0, len(lines)):
                        if lines[i] == '[Bioinformatics]':
                            # when Bioinformatics section is found, start
                            # copying lines to the list buffer, but don't
                            # copy the header itself.
                            sentinel = True
                        elif lines[i].startswith('['):
                            # when the next header is found, stop copying
                            # lines to the list buffer. Don't include this
                            # header line, either.
                            sentinel = False
                        elif sentinel == True:
                            # this should be a line in between
                            # [Bioinformatics] and the next section. Copy it
                            # to the buffer.
                            metadata.append(lines[i])
                        # if the end of the file is reached before the next
                        # header is found, that means there were no more
                        # headers and that's okay too.

                    # remove duplicate lines (this appears to be an issue in
                    # the original bash scripts.)
                    l = list(set(metadata))
                    l.sort()

                    # write the sanitized data out to legacy file.
                    # TODO: Are processing these one at a time, or otherwise can allow overwriting the file
                    #  for each CSV found?
                    with open(self.run_config_file_path, 'w') as f2:
                        for line in metadata:
                            # end lines w/proper UNIX-style newline, unless
                            # Windows '\r\n' is preferred.
                            f2.write("%s\n" % line)

    def run(self):
        if not os.path.exists(self.run_config_file_path):
            # run_param doesn't run anything, it actually generates the run_config.txt file. :)
            self._run_param()

        fastq_output = "%s/Data/Fastq" % self.seq_dir

        '''
        sample contents of run_config.txt files:
        
        Project,ForwardAdapter,ReverseAdapter,PolyGTrimming,HumanFiltering,QiitaID
        apple_test_1,NA,NA,FALSE,FALSE,10556
        apple_test_1,NA,NA,FALSE,FALSE,10567
        apple_test_1,NA,NA,FALSE,FALSE,10568
        apple_test_1,NA,NA,FALSE,FALSE,10569
        apple_test_1,NA,NA,FALSE,FALSE,10570
        apple_test_1,NA,NA,FALSE,FALSE,10571

        and

        project_name,ForwardAdapter,ReverseAdapter,PolyGTrimming,HumanFiltering,QiitaID
        iseq_44444,GATCGGAAGAGCACACGTCTGAACTCCAGTCAC,GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT,TRUE,FALSE,44444
        '''





'''
  while IFS=" " read project adapter_a adapter_A a_trim h_filter qiita_proj; do


    proj_check=$(echo $project | awk -F_ '{print $NF}')

    ### remove any trailing newlines (^M) from variable
    qiita_proj=$(echo -n $qiita_proj)
    qiita_proj=$(echo "$qiita_proj"|tr -d '\r')
    pushd ${fastq_output}/${project}
      input_count=$(ls ${fastq_output}/${project}/*R1*.fastq.gz | wc -l)

      if [[ "${input_count}" -gt "2000" ]]; then
        split_count=16
      elif [[ "${input_count}" -le "2000" && "${input_count}" -gt "1000" ]]; then
        split_count=10
      elif [[ "${input_count}" -le "1000" && "${input_count}" -gt "500" ]]; then
        split_count=4
      else
        split_count=1
      fi

      if [[ -z $trim_file ]]; then
      	trim_file=split_file_
      fi

      job_count=$(($split_count + 1))
      file_base=$(basename $seqdir)

      if [[ -e ${trim_file}* ]]; then
        rm -fv ${trim_file}*
      fi

      find . -maxdepth 1 -name "*.fastq.gz" -type f | grep "_R1_" | cut -f2 -d"/" > ${file_base}_file_list.txt
      #ls | grep *.fastq.gz > ${file_base}_file_list.txt
      line_count=$(cat ${file_base}_file_list.txt | wc -l)
      line_per_split=$(( $(( $line_count + $split_count -1 ))/$split_count ))

      split -l $line_per_split -d ${file_base}_file_list.txt split_file_; rename ${trim_file}0 ${trim_file} ${trim_file}*
      if [[ ${line_count} -eq "0" ]]; then
        :
      fi

    NOTE:  Create ${seqdir}/Data/Fastq/${project}/${file_base}_qsub.sh


    pbs_job_id=$(sbatch --qos=seq_proc --parsable --array=0-$(($split_count - 1)) ${seqdir}/Data/Fastq/${project}/${file_base}_qsub.sh)

    ### this should never fail as it only is executed while there is another project to process
    ### first time through
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

}


'''
