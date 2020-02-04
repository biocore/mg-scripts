#!/bin/bash

### usage:
### sh bowtie_trime_qsub.sh <directory to filter> <split file name prefix> <number of files to split into>

set -x

if [[ $# != 3 ]]; then
	echo not enough arguments
	echo Usage: sh bowtie_trim_qsub.sh source directory split_file_prefix num_split_files
	exit 1
fi

if [[ ! -d $1 ]]; then
	echo source directory $1 does not exist
	exit 2
fi

if [[ $3 > 8 ]]; then
	echo script not built to handle $3 split files
	echo please select split file < 9
	exit 3
fi


dir=$1
trim_file=$2
split_count=$3 ### must be less than 9
job_count=$(($split_count + 1))

file_base=$(basename $dir)

pushd $dir

find . -maxdepth 1 -name "*.fastq.gz" -type f | grep "_R1_" | cut -f2 -d"/" > ${file_base}_file_list.txt
line_count=$(cat ${file_base}_file_list.txt | wc -l)

split -l $(( $line_count / 5)) -d ${file_base}_file_list.txt split_file_

if [[ $? == 0 ]]; then
	rm ${file_base}_file_list.txt
fi


#qsub -v dir="$dir",trim_file="$trim_file" -t 0-$(($split_count + 1)) ~/filter_job_parallel.sh
pbs_job_id=$(qsub -v dir="$dir",trim_file="$trim_file" -t 0-$(($split_count + 1)) ~/filter_job_parallel.sh)

qsub -v dir="$dir" -lnodes=1:ppn=16 -lpmem=1gb ~/fastqc_parallel.sh -W depend=afterok:$pbs_job_id
 

