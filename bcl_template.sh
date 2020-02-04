#!/bin/bash

### #PBS -N $expname 
###bcl2fastq
#PBS -l nodes=1:ppn=14
#PBS -V
#PBS -l walltime=48:00:00
#PBS -l pmem=10G
#PBS -o ${job_o_out} 
###/sequencing/sequencer_mount
#PBS -e ${job_e_out}

#PBS -m abe
#PBS -M ${user_email}

dataprefix=${final_data_location}
bcl2fastq="~/bcl2fastq_2_20/bin/bcl2fastq"

module load gcc_4.9.1

#dataloc=$1
#outputdir=$2
if [[ ! -z ${base_mask} ]]; then
	use_base_mask="${base_mask}"
fi

cd ${seqdir}

#if [ "x$PBS_O_WORKDIR" != "x" ]; then
#        cd $PBS_O_WORKDIR
#fi

export PATH=$PATH:/usr/local/bin

#sleep 300

cmd="$bcl2fastq \
--sample-sheet ${csvfile} \
--minimum-trimmed-read-length 1 \
--mask-short-adapter-reads 1 \
-R . \
-o ${outputdir} \
--loading-threads 4 \
--processing-threads 8 \
--writing-threads 2 \
--ignore-missing-bcls \
--create-fastq-for-index-reads "
#--barcode-mismatches 0 \
#--ignore-missing-bcls" 
#--tiles s_[12] --ignore-missing-bcls"  
#--barcode-mismatches 0"
#--ignore-missing-bcls"

date
echo "Executing: $cmd"
eval $cmd
date

