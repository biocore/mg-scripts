#!/bin/bash
set -x

source ~/.seqvars

###source /home/jede9131/miniconda3/bin/activate test_env_2
module load fastp_0.20.1 samtools_1.12 minimap2_2.18

while getopts ":d:D:S:s:p:c:B:o:O:a:A:q:g:G:C:m:f:" opt; do
    case ${opt} in
        x ) #path to human-phix-db.mmi database file
            db=$OPTARG
            ;;
        d ) #directory to process
            dir=$OPTARG
            ;;
        D ) ##location for fastq files
            file_dir=$OPTARG
            ;;
        S ) #trim_file prefix
            trim_file=$OPTARG
            ;;
        s ) # number of files to split input
            split_count=$OPTARG
            ;;
        p ) #specific project to process
            project=$OPTARG
            ;;
        c ) # cores to use
            NPROCS=$OPTARG
            ;;
        B ) ### file with files to reprocess
            bad_file=$OPTARG
            ;;
        o ) #optional outputdir for debug
            ### .seqvars outputdir + run_identifier
            output_dir=$OPTARG
            ;;
        O ) ### run identifer
            base_seq_dir=$OPTARG
            ;;
        a ) #adapter 1
            adapter_a=$OPTARG
            ;;
        A ) #adapter 2
            adapter_A=$OPTARG
            ;;
        q ) #qiita projects
            qiita_proj=$OPTARG
            ;;
        g ) #adapter trimming
            a_trim=$OPTARG
            ;;
        G ) #human filter
            h_filter=$OPTARG
            ;;
        C ) #Chemistry
            chemistry=$OPTARG
            ;;
        m ) ### copy files?
            copy_file=$OPTARG
            ;;
        f ) ### final output Directory
            final_output_dir=$OPTARG
            ;;
        \? ) echo "Usage: cmd [-d <directory>] [-t <split file prefix>] [-s <num split files>] [-p <specific project>] [-c cores =16]"
            ;;
        : ) "$OPTARG requires argument" 1>&2
            ;;
    esac
done
shift $((OPTIND -1))

if [[ -z $dir ]]; then
    echo Directory is required
    exit 1
fi
### set a couple defaults
if [[ -z $NPROCS ]]; then
    NPROCS=16
fi
if [[ -z $trim_file ]]; then
    trim_file=split_file_
fi
if [[ -z $output_dir ]]; then
    output_dir=filtered_sequences_test_2
fi
if [[ -z $a_trim ]]; then
    a_trim=false
fi
if [[ -z $h_filter ]]; then
    h_filter=false
fi
if [[ -z $qiita_proj ]]; then
    qiita_proj="NA"
fi
if [[ -z $project ]]; then
    grep_var="egrep -v 'Stats|Reports|moved|sample|xtract'"
elif [[ $project ]]; then
    grep_var="grep $project"
fi
if [[ -z $chemistry ]]; then
    chemistry="Default"
fi
### copy files to qiita by default if project defined.
if [[ -z $copy_file ]]; then
    copy_file="TRUE"
fi

data_prefix=/var/www/html/seq_fastqc
#bowtie=$(which bowtie2)
samtools=$(which samtools)
#bedtools=$(which bedtools)
fastp=$(which fastp)
minimap2=$(which minimap2)
gui_output=gui_files
indiv_output=run_files


suffix=R*.fastq*

if [[ -z $NPROCS ]]; then
    NPROCS=16
fi

### parallel gzip??
if [[ `which pigz` ]]; then
    _zip="pigz -p $NPROCS -f"
else
    _zip="gzip -f"
fi

#filter_db="/databases/bowtie/Human_phiX174/Human_phix174"
#filter_db="/databases/bowtie/phiX/phiX"
tar="/bin/tar"

fastp_qc_output=${dir}/${project}/fastp_qc

if [[ -e ${rundir}/run_config.sh ]]; then
    source ${rundir}/run_config.sh
    ### parse runconfig file
fi

### default settings
#if [[ -z $adapter_a ]]; then
#  adapter_a="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
#fi
#if [[ -z $adapter_A ]]; then
#  adapter_A="GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT"
#fi
#if [[ -z $adapter_polyg_a ]]; then
#	adapter_polyg_a="GGGGGGGGGG"
#fi
#if [[ -z $adapter_polyg_A ]]; then
#	adapter_polyg_A="GGGGGGGGGG"
#fi
if [[ -z $filter_db ]]; then
    #  filter_db="/databases/bowtie/Human_phiX174/Human_phix174"
    filter_db="/databases/bowtie/Human/Human"
fi


dir=${dir}/${file_dir}
pushd $dir/$project

final_output=$output_dir/${base_seq_dir}

if [[ ! -d ${final_output} ]]; then
    mkdir ${final_output}

fi

if [[ ! -d ${final_output}/${project} ]]; then
    mkdir ${final_output}/${project}
fi

### move sample sheets into place
### pull everything from original location to top level output
cp ${dir}/*.csv $output_dir

dir_arr=("html" "json")
for stat_dir in "${dir_arr[@]}"; do
    if [[ ! -d ${final_output}/${project}/$stat_dir ]]; then
        mkdir ${final_output}/${project}/$stat_dir
        mkdir ${dir}/${project}/${stat_dir}
    fi
done
### moved to individual tri/mmed/filtered directory
#if [[ -d ${final_output}/${location}/filtered_sequences ]]; then
#	mkdir ${final_output}/${location}/filtered_sequences
#fi
dir_arr=("trim_logs" "zero_files")
for location in "${dir_arr[@]}"; do
    if [[ ! -d ${final_output}/${project}/${location} ]]; then
        mkdir ${final_output}/${project}/${location}
    fi
done

if [[ ! -d ${seqdir} ]]; then
    mkdir -p $fastp_qc_output
    mkdir ${fastp_qc_output}/fastp_logs
    mkdir ${fastp_qc_output}/fastp_fastqc
fi


if [[ ! -e $trim_file ]]; then
    find . -maxdepth 1 -name "*.fastq.gz" -type f | grep "_R1_" | cut -f2 -d"/" > ${trim_file}
fi

### convert true/false/yes/no to uppercase
a_trim=$(echo ${a_trim^^})
h_filter=$(echo ${h_filter^^})

### amplicon MiSeq sequencing
if [[ ($a_trim =~ ^(FALSE|NO)$ ) && ($h_filter =~ ^(FALSE|NO)$ ) ]]; then ### && ( ! -z $qiita_proj) ]]; then
    ### just copy data to qiita study.  16s data.
    ###if [[ $qiita_proj ]]; then
    if [[ $chemistry == "Amplicon" ]]; then
        mkdir ${final_output_dir}/${project}/amplicon
        ### copy all files to final location
        cp ${dir}/${project}/*.fastq.gz ${final_output_dir}/${project}/amplicon
        ###:
    fi
    if [[ $qiita_proj ]]; then
        filename_index_1=$(echo "$filename1" | sed -e 's/_R1_/_I1_/g')
        filename_index_2=$(echo "$filename1" | sed -e 's/_R2_/_I2_/g')
        ### modify this to pull fastq.gz files from original $seqdir location.  this is where they
        ### would have been written from bcl2fastq.  no trimming or filtering would have been done
        ### so no trimmed/filtered final_output/project locations.
        ###if [[  ! $copy_file ]]; then
        ###rsync -avp --progress --rsync-path="pfexec rsync" --chown qiita:qiita ${dir}/${project}/*.fastq.gz jdereus@10.210.38.14:/z2pool-1/data/qiita_data/uploads/${qiita_proj}/
        if [[ $qiita_proj != "NA" &&  $copy_file != "FALSE" ]]; then
            pushd ${dir}/${project}
            sudo -u qiita rsync -avp --progress *.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}
        fi
        #if [[ $(stat -c '%U' /qmounts/qiita_data/uploads/${qiita_proj} ) != 'qiita' ]]; then
        #	sudo chown -v 5500:5500 /qmounts/qiita_data/uploads/${qiita_proj}
        #fi
        popd

        ###fi
        ###			rsync -avp --progress --rsync-path="pfexec rsync" --chown qiita:qiita ${final_output}/${project}/trimmed_sequences/${transfer_3}.fastq.gz ${final_output}/${project}/trimmed_sequences/${transfer_4}.fastq.gz jdereus@10.210.38.14:/z2pool-1/data/qiita_data/uploads/${qiita_proj}/
    fi
fi

###possible amplicon
if [[ ($a_trim =~ ^(TRUE|YES)$ ) && ($h_filter =~ ^(FALSE|NO)$ ) ]]; then
    if [[ ! -d ${final_output}/${project}/trimmed_sequences ]]; then
        mkdir -p ${final_output}/${project}/trimmed_sequences
    fi
    echo "here i am inside true/false"
    for file in `cat ${dir}/${project}/${trim_file} | grep "_R1_"`;
    do
        parent_dir=$(dirname $file)
        project_dir=$(echo $parent_dir | cut -f 2 -d"/")
        ### actual filename minus the preceding path
        filename1=$(basename "$file") ### .fastq.gz)
        ### strip fastq.gz for better file naming regarding output
        filename1_short=$(basename "$filename1" .fastq.gz)
        ### replace R1 with R2 to be able to pass both forward and reverse reads to bowtie
        filename2=$(echo "$filename1" | sed -e 's/_R1_00/_R2_00/g')
        filename2_short=$(basename "$filename2" .fastq.gz)
        filename_index1=$(echo "$filename1" | sed -e 's/_R1_/_I1_/g')
        filename_index2=$(echo "$filename2" | sed -e 's/_R2_/_I2_/g')

        ### changed output of json/html from ${dir} to ${final_output} 10/6/20
        #if [[ $adapter_a == 'NA' || $adapter_A == 'NA' ]]; then
        if [[ -z $adapter_a || -z $adapter_A || $adapter_a == "NA" || $adapter_A == "NA" ]]; then
            $fastp -l 100 -i $filename1 -I $filename2 -w $NPROCS -j ${final_output}/${project}/json/${filename1_short}.json -h ${final_output}/${project}/html/${filename1_short}.html -o $final_output/${project}/trimmed_sequences/${filename1_short}.fastp.fastq.gz -O $final_output/${project}/trimmed_sequences/${filename2_short}.fastp.fastq.gz -R ${filename1_short}_report
        else
            $fastp --adapter_sequence ${adapter_a} --adapter_sequence_r2 ${adapter_A} -l 100 -i $filename1 -I $filename2 -w $NPROCS -j ${final_output}/${project}/json/${filename1_short}.json -h ${final_output}/${project}/html/${filename1_short}.html -o $final_output/${project}/trimmed_sequences/${filename1_short}.fastp.fastq.gz -O $final_output/${project}/trimmed_sequences/${filename2_short}.fastp.fastq.gz -R ${filename1_short}_report ### -h ${filename1_short}.html
        fi

        #pushd $final_output/${project}/trimmed_sequences/

        #$_zip ${filename1_short}.fastp.fastq
        #$_zip ${filename2_short}.fastp.fastq
        #popd

        size1=$(stat -c%s $final_output/${project}/trimmed_sequences/${filename1_short}.fastp.fastq.gz)
        size2=$(stat -c%s $final_output/${project}/trimmed_sequences/${filename2_short}.fastp.fastq.gz)

        ### amplicon/16s data. include index files
        if [[ $chemistry == "Amplicon" ]]; then
            cp ${filename_index1}* $final_output/${project}/trimmed_sequences
            cp ${filename_index2}* $final_output/${project}/trimmed_sequences
        fi

        if [[ $size1 -le 500 || $size2 -le 500 ]]; then
            echo $filename1_short is very small
            mv $final_output/${project}/trimmed_sequences/${filename1_short}* $final_output/${project}/trimmed_sequences/${filename2_short}* ${final_output}/${project}/zero_files
            echo $final_output/${project}/trimmed_sequences${filename1_short}.fastp.fastq.gz >> $final_output/${project}/empty_file_list.txt
            echo $size1 >> $final_output/${project}/empty_file_list.txt
            echo $final_output/${project}/trimmed_sequences/${filename2_short}.fastp.fastq.gz >> $final_output/${project}/empty_file_list.txt
            echo $size2 >> $final_output/${project}/empty_file_list.txt
        else
            :
            ### transfer each file upon completion?
            #rsync -ap -e "ssh -i ${HOME}/.ssh/qiita_rsa" ${final_output}/${project}/${filename1_short}.fastp.fastq.gz ${final_output}/${project}/${filename2_short}.fastp.fastq.gz qiita@10.210.38.7:/projects/qiita_data/uploads/${qiita_proj}/
        fi

    done
    ### done.  move files to qiita study if possible
    echo qiita project == $qiita_proj

    #if [[ ! $copy_file ]]; then
    ### qiita project is defined
    #if [[ $qiita_proj != "NA" ]]; then
    for file_transfer in `cat ${dir}/${project}/${trim_file}`; do
        transfer_1=$(echo "$file_transfer" | cut -f1 -d".")
        transfer_2=$(echo "$transfer_1" | sed -e 's/_R1_00/_R2_00/g')
        transfer_3=$(echo "$transfer_1" | sed -e 's/_R1_00/_I1_00/g')
        transfer_4=$(echo "$transfer_1" | sed -e 's/_R1_00/_I2_00/g')
        ### uncomment the following for panfs locations
        sleep 1
        if [[ $chemistry == "Amplicon" ]]; then
            ### copy all files
            rsync -avp --progress ${final_output}/${project}/*.fastp.fastq.gz ${final_output_dir}/${base_seq_dir}/${project}
        else
            rsync -avp --progress ${final_output}/${project}/${transfer_1}.fastp.fastq.gz ${final_output}/${project}/${transfer_2}.fastp.fastq.gz ${final_output_dir}/${base_seq_dir}/${project}
        fi
        #sudo -u qiita rsync -avp --progress --chown qiita:qiita ${final_output}/${project}/trimmed_sequences/${transfer_1}.fastp.fastq.gz ${final_output}/${project}/trimmed_sequences/${transfer_2}.fastp.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/
        ###rsync -avp --progress --rsync-path="pfexec rsync" --chown qiita:qiita ${final_output}/${project}/trimmed_sequences/${transfer_1}.fastp.fastq.gz ${final_output}/${project}/trimmed_sequences/${transfer_2}.fastp.fastq.gz jdereus@10.210.38.14:/z2pool-1/data/qiita_data/uploads/${qiita_proj}/
        ### if amplicon copy I1/I2 files to project directory

        if [[ $qiita_proj != "NA" ]]; then
            if [[ $chemistry == "Amplicon" && $copy_file == "TRUE" ]]; then
                sudo -u qiita rsync -avp --progress --chown 5500:5500 ${final_output}/${project}/trimmed_sequences/*_R[12]*.fastq.gz ${final_output}/${project}/trimmed_sequences/*_I1*.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/
                ###${final_output}/${project}/trimmed_sequences/${transfer_3}.fastq.gz ${final_output}/${project}/trimmed_sequences/${transfer_4}.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/
                ###sudo -u qiita rsync -avp --progress --chown 5500:5500 ${final_output}/${project}/trimmed_sequences/${transfer_1}.fastq.gz ${final_output}/${project}/trimmed_sequences/${transfer_2}.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/
            elif [[ $copy_file == "TRUE" ]]; then
                sudo -u qiita rsync -avp --chown 5500:5500 ${final_output}/${project}/trimmed_sequences/${transfer_1}.fastp.fastq.gz ${final_output}/${project}/trimmed_sequences/${transfer_2}.fastp.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/
            fi
        fi
        ### if [[ $(stat -c '%U' /qmounts/qiita_data/uploads/${qiita_proj} ) != 'qiita' ]]; then
        ###sudo chown -v 5500:5500 /qmounts/qiita_data/uploads/${qiita_proj}
        ###ssh jdereus@10.210.38.14 "pfexec /usr/gnu/bin/chown -v qiita:qiita /z2pool-1/data/qiita_data/uploads/${qiita_proj}/"
        ##fi
    done
    sleep 1
    rsync -avp --progress ${final_output}/${project}/ ${final_output_dir}/${base_seq_dir}
    #fi
    #fi

elif [[ ($a_trim =~ ^(TRUE|YES)$) && ($h_filter =~ ^(TRUE|YES)$) ]]; then
    echo "here i am inside true/true"
    if [[ ! -d ${final_output}/${project}/filtered_sequences ]]; then
        mkdir ${final_output}/${project}/filtered_sequences
    fi

    for file in `cat ${dir}/${project}/${trim_file} | grep "_R1_"`;
    do
        parent_dir=$(dirname $file)
        project_dir=$(echo $parent_dir | cut -f 2 -d"/")
        ### actual filename minus the preceding path
        filename1=$(basename "$file") ### .fastq.gz)
        ### strip fastq.gz for better file naming regarding output
        filename1_short=$(basename "$filename1" .fastq.gz)
        ### replace R1 with R2 to be able to pass both forward and reverse reads to bowtie
        filename2=$(echo "$filename1" | sed -e 's/_R1_00/_R2_00/g')
        filename2_short=$(basename "$filename2" .fastq.gz)

        if [[ $adapter_a == "NA" || $adapter_A == "NA" ]]; then
            $fastp -l 100 -i $filename1 -I $filename2 -w $NPROCS --stdout -j ${dir}/${project}/json/${filename1_short}.json -h ${dir}/${project}/html/${filename1_short}.html | $minimap2 -ax sr -t $NPROCS $db - -a | $samtools fastq -@ $NPROCS -f 12 -F 256 -1 $final_output/${project}/filtered_sequences/${filename1_short}.trimmed.fastq.gz -2 $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz
        else
            $fastp --adapter_sequence ${adapter_a} --adapter_sequence_r2 ${adapter_A} -l 100 -i $filename1 -I $filename2 -w $NPROCS --stdout -j ${dir}/${project}/json/${filename1_short}.json -h ${dir}/${project}/html/${filename1_short}.html | $minimap2 -ax sr -t $NPROCS $db - -a | $samtools fastq -@ $NPROCS -f 12 -F 256 -1 $final_output/${project}/filtered_sequences/${filename1_short}.trimmed.fastq.gz -2 $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz
        fi
        size1=$(stat -c%s $final_output/${project}/filtered_sequences/${filename1_short}.trimmed.fastq.gz)
        size2=$(stat -c%s $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz)

        if [[ $size1 -le 500 || $size2 -le 500 ]]; then
            mv $final_output/${project}/filtered_sequences/${filename1_short}* $final_output/${project}/filtered_sequences/${filename2_short}* ${final_output}/${project}/zero_files
            echo $final_output/${project}/filtered_sequences${filename1_short}.trimmed.fastq.gz >> $final_output/${project}/empty_file_list.txt
            echo $size1 >> $final_output/filtered_sequences/${project}/empty_file_list.txt
            echo $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz >> $final_output/${project}/empty_file_list.txt
            echo $size2 >> $final_output/filtered_sequences/${project}/empty_file_list.txt
        else
            :
            ### transfer each file upon completion?
            #rsync -ap -e "ssh -i ${HOME}/.ssh/qiita_rsa" ${final_output}/${project}/${filename1_short}.fastp.fastq.gz ${final_output}/${project}/${filename2_short}.fastp.fastq.gz qiita@10.210.38.7:/projects/qiita_data/uploads/${qiita_proj}/
        fi

    done

    ### always executes
    if [[ $copy_file ]]; then
        for file_transfer in `cat ${dir}/${project}/${trim_file}`; do
            transfer_1=$(echo "$file_transfer" | cut -f1 -d".")
            transfer_2=$(echo "$transfer_1" | sed -e 's/_R1_00/_R2_00/g')
            rsync -avp --progress ${final_output}/${project}/filtered_sequences/${transfer_1}.trimmed.fastq.gz ${final_output}/${project}/filtered_sequences/${transfer_2}.trimmed.fastq.gz ${final_output_dir}/${base_seq_dir}/${project}/filtered_sequences/

            if [[ $qiita_proj != "NA" && $copy_file == "TRUE" ]]; then
                sudo -u qiita rsync -avp --chown 5500:5500 ${final_output}/${project}/filtered_sequences/${transfer_1}.trimmed.fastq.gz ${final_output}/${project}/filtered_sequences/${transfer_2}.trimmed.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/
            fi ###rsync -avp --rsync-path="pfexec rsync" --chown qiita:qiita ${final_output}/${project}/filtered_sequences/${transfer_1}.trimmed.fastq.gz ${final_output}/${project}/filtered_sequences/${transfer_2}.trimmed.fastq.gz jdereus@10.210.38.14:/z2pool-1/data/qiita_data/uploads/${qiita_proj}/
        done
        rsync -avp --progress ${final_output} ${final_output_dir}/
        #			if [[ $(stat -c '%U' /qmounts/qiita_data/uploads/${qiita_proj} ) != 'qiita' ]]; then
        #				sudo chown -v 5500:5500 /qmounts/qiita_data/uploads/${qiita_proj}/
        #			fi
        ###ssh jdereus@10.210.38.14 "pfexec /usr/gnu/bin/chown -v qiita:qiita /z2pool-1/data/qiita_data/uploads/${qiita_proj}/"
        ###fi
        ###fi

    fi
fi


if [[ $h_filter =~ ^(FALSE|NO)$ ]]; then
    ### we are done.  copy adapter trimmed to qiita_project
    : ###rsync -ap --rsync-path="pfexec rsync" ${fastp_qc_output}/*R[12].fastq.gz jdereus@10.210.38.14:/z2pool-1/data/qiita_data/uploads/${qiita_proj}/
    if [[ $? -ne "0" ]]; then
        echo "transfer failed ${dir}" | mailx -s "failed transfer" jdereus@ucsd.edu
        exit -1
    else
        exit 0
    fi
else
    echo "moving on to human filtering"
    :
fi
#  echo "-name "*$suffix" -type f -exec fastqc {} -t ${NPROCS} -o $fastp_qc_output \;"
#	find $fastp_qc_output -name "*$suffix" -maxdepth 1 -type f -exec fastqc {} -t ${NPROCS} -o $fastp_qc_output/fastp_fastqc \;
#	find $fastq_dir -name "*$suffix" -maxdepth 1 -type f -exec fastqc {} -t ${NPROCS} -o $output_dir \;

#  cd $fastp_qc_output/fastp_fastqc
#  multiqc ./

### pop back up to $dir
#popd

###pushd ${fastp_qc_output}


#####################
bowtie shouldn't be called. minimap should be doing the human filtering
conditionals should make it so we exit with 0 or -1 and don't move onto using bowtie
bowtie + atropos = wrong. minimap + fastp = good.
send an email to jeff to make sure that bowtie and atropos were never called. if so we need to know when and why.
#####################

#for file in `find ${fastp_qc_output} -maxdepth 1 -type f -name "*.fastq" | grep _R1_`;
#for file in `cat ${fastp_qc_output}/${trim_file}`;
if [[ $h_filter == FALSE-ISH ]]; then
    exit 0
    ### explicit check for now
elif [[ $h_filter == TRUE-ISH ]]; then
    for file in `cat ${dir}/${project}/${trim_file}`;
    do
        final_output=${dir}/${project}/filtered_sequences
        if [[ ! -d ${final_output} ]]; then
            mkdir ${final_output}
            mkdir ${final_output}/trim_logs
            mkdir ${final_output}/zero_files
        fi

        if [[ ! -d ${fastp_qc_output}/processed_qc ]]; then
            mkdir ${fastp_qc_output}/processed_qc
        fi


        parent_dir=$(dirname $file)
        project_dir=$(echo $parent_dir | cut -f 2 -d"/")
        ### actual filename minus the preceding path
        filename1=$(basename "$file") ### .gz) ### .fastq.gz)
        ### strip fastq.gz for better file naming regarding output
        filename1_short=$(basename "$filename1" .fastq.gz)
        filename1_trim=$(echo "filename1_short" | sed -e 's/atropos//g')
        ### replace R1 with R2 to be able to pass both forward and reverse reads to bowtie
        filename2=$(echo "$filename1" | sed -e 's/_R1_/_R2_/g')
        filename2_short=$(basename "$filename2" .fastq.gz)
        filename2_trim=$(echo "filename2_short" | sed -e 's/atropos//g')

        $bowtie -p $NPROCS -x ${filter_db} --very-sensitive -1 ${fastp_qc_output}/$filename1 -2 ${fastp_qc_output}/$filename2 2> ${final_output}/trim_logs/${filename1_short}.log | $samtools view -f 12 -F 256 | $samtools sort -@ 16 -n | $samtools view -bS | $bedtools bamtofastq -i - -fq $final_output/${filename1_short}.trimmed.fastq -fq2 $final_output/${filename2_short}.trimmed.fastq &> $final_output/trim_logs/${filename1_short}.log
        #        gzip ${final_output}/${filename1_short}.filtered.fastq.gz ${final_output}/${filename1_short}.trimmed.fastq
        #        gzip ${final_output}/${filename2_short}.filtered.fastq.gz ${final_output}/${filename2_short}.trimmed.fastq
        #		$bowtie -p 16 -x ${filter_db} --very-sensitive -1 ${fastp_qc_output}/$filename1 -2 ${fastp_qc_output}/$filename2 | $samtools view -f 12 -F 256 | $samtools sort -@ 16 -n | $samtools view -bS | $bedtools bamtofastq -i - -fq $final_output/${filename1_short}.fastq -fq2 $final_output/${filename2_short}.fastq &> $final_output/trim_logs/${filename1_short}.log
        $_zip ${final_output}/${filename1_short}.trimmed.fastq
        $_zip ${final_output}/${filename2_short}.trimmed.fastq

        if [[ $(gzip -l ${final_output}/${finalname1_short}.trimmed.fastq.gz | awk 'NR==2 {exit($2!=0)}') ]]; then
            echo ${final_output}/${finalname1_short}.trimmed.fastq.gz >> zero_files.txt
            rm ${final_output}/${finalname1_short}.trimmed.fastq.gz ${final_output}/${finalname2_short}.trimmed.fastq.gz
            #		mv ${final_output}/${finalname1_short}.trimmed.fastq.gz ${final_output}/${finalname2_short}.trimmed.fastq.gz ${final_output}/zero_files
        else
            :
        fi

        for file in `cat ${dir}/${project}/${trim_file} | cut -f 1,2,3 -d "_"`;
        do
            if [[ $qiita_proj != "NA" ]]; then
                #file=$(file | cut -f 1,2,3 -d "_")
                ### copy files to qiita study
                echo copying files ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz
                #rsync -ap --rsync-path="pfexec rsync" ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz jdereus@10.210.38.14:/z2pool-1/data/qiita_data/uploads/${qiita_proj}/
            else
                ### samplesheet has a qiita study ID listed
                echo syncing file ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz
                #rsync -ap ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz $output_dir/${project}/
            fi
        done

    done
fi

popd
