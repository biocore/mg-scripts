

class FastQCParallelBCLConvert:
    def __init__(self, dir, output_dir, nprocs, proj_dir, lab_name, chemistry):
        # source ~/.seqvars if needed.

        if chemistry.lower() == 'amplicon':
            self.suffix = "[RI]*.fastq*"
        else:
            self.suffix = "R*.fastq*"

        data_prefix = "/projects/fastqc"
        # source /home/jede9131/miniconda3/bin/activate test_env_2

    def run_param(self):
        # for csvfile in $(ls ${dir}/*.csv | grep -v "sav.csv"); do
        #     sed -n '/\[Project/,/Contact/ p' ${dir}/${csvfile} | egrep -v 'Project|Contact' > ${dir}/run_config.txt
        #     sed -i 's/,/ /g' ${dir}/run_config.txt
        #     sed -i '/^\s*$/d' ${dir}/run_config.txt
        # done
        pass

    def do_work(self):
        pass
        # pushd $dir
        # final_output=$output_dir
        #
        # ### removed this converstion 10/6/20
        # #output_dir=$(basename $dir)
        # fastq_raw=${dir}/Data/Fastq
        # atropos_qc_output=${dir}/${proj_dir}/atropos_qc/processed_qc
        # #fastq_trimmed=${dir}/Data/Fastq ###/${proj_dir}/filtered_sequences
        #
        # fastq_trimmed=${output_dir}/$(basename ${dir})/${proj_dir}
        #
        # parallel=`which parallel`
        # tree=`which tree`
        #
        # ### running on single directory out of pipeline
        # if [[ ! -f ${dir}/run_config.txt ]]; then
        #   run_param
        # fi
        #
        # while read line; do
        #   IFS=' ' proj_arr+=( $(echo ${line} | awk '{print $1}') )
        # done < ${dir}/run_config.txt
        #
        # fastqc_output_dir=${labname}/$(basename `pwd`)
        # if [[ "mount | grep $data_prefix" ]]; then
        #   : #statements
        # else
        #   exit 666 ### no network directory
        # fi
        #
        # initial="true"
        # if [[ $initial == "true" ]]; then
        #   for project in $proj_dir; do ###${proj_arr[@]}; do
        #
        #     if [[ $(find ${fastq_trimmed}/ -maxdepth 2 -name "trimmed_sequences" -type d) ]]; then
        #       filter_type=trimmed_sequences
        #     filter_type_sub=trimmed
        #     elif [[ $(find ${fastq_trimmed}/ -maxdepth 2 -name "filtered_sequences" -type d) ]]; then
        #       filter_type=filtered_sequences
        #     filter_type_sub=filtered
        #     elif [[ $(find ${fastq_trimmed}/ -maxdepth 2 -name "amplicon" -type d) ]]; then
        #       filter_type=amplicon
        #     filter_type_sub=amplicon
        #     fi
        #
        #     mkdir -p ${data_prefix}/${fastqc_output_dir}/${project}/{bclconvert,$filter_type}
        #   find ${fastq_raw}/${project} -maxdepth 1 -name  "*$suffix" -type f | ${parallel} -j $NPROCS fastqc --noextract -t ${NPROCS} {} -o ${data_prefix}/${fastqc_output_dir}/${project}/bclconvert/ \;
        #
        #
        #   pushd ${data_prefix}/${fastqc_output_dir}/${project}/bclconvert/
        #   popd
        # popd
        #
        #
        # pushd ${fastq_trimmed}
        #   find ${fastq_trimmed}/${filter_type}/  -maxdepth 1 -name "*$suffix" -type f | ${parallel} -j ${NPROCS} fastqc --noextract -t ${NPROCS} {}  -o ${data_prefix}/${fastqc_output_dir}/${project}/${filter_type} \;
        #   pushd ${data_prefix}/${fastqc_output_dir}/${project}/${filter_type}
        #   popd
        # popd
        #
        #
        # if [[ ! -d ${data_prefix}/${fastqc_output_dir}/${proj_dir}/multiqc ]]; then
        #   mkdir -p ${data_prefix}/${fastqc_output_dir}/${proj_dir}/multiqc
        # fi
        # pushd ${dir}
        # ### remove Stats directory.  Stats is bcl2fastq specific.  bclconvert writes to Reports.
        #   multiqc -c ~/multiqc-bclconvert-config.yaml --fullnames --force ${data_prefix}/${fastqc_output_dir}/Reports ${fastq_raw}/${proj_dir}/json ${fastq_trimmed}/json/ ${data_prefix}/${fastqc_output_dir}/${proj_dir}/${filter_type} ${data_prefix}/${fastqc_output_dir}/${proj_dir}/bclconvert -o ${data_prefix}/${fastqc_output_dir}/${proj_dir}/multiqc --interactive
        # popd
        # done
        #
        # pushd ${data_prefix}
        #   $tree -H '.' -d -t -L 1 --noreport --charset utf-8 > index.html
        # popd
        #
        # pushd ${data_prefix}/${fastqc_output_dir}/
        #   $tree -H '.' -d -t -L 1 --noreport --charset utf-8 > index.html
        # popd
        # pushd ${data_prefix}/${labname}
        #   $tree -H '.' -d -t -L 1 --noreport --charset utf-8 > index.html
        # popd
        #
        # pushd ${data_prefix}/${fastqc_output_dir}/${proj_dir}/multiqc
        # popd
        #
        # fi
        # echo result==$?
        # if [[ $? == 0 ]]; then
        #   ( printf '%s\n\n' "fastqc complete for project ${project}"
        #     printf '%s\n\n' "fastqc results located at http://kl-fastqc.ucsd.edu/fastqc/${fastqc_output_dir}/${project}"
        #   ) | mailx -v -s "fastqc complete ${project}" -r jdereus@ucsd.edu $fastqc_mail
        # fi









