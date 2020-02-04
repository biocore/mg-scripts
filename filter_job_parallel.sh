#!/bin/bash

### #PBS -N $expname 
#PBS -l nodes=1:ppn=16
#PBS -V
#PBS -l walltime=384:00:00 
#PBS -l pmem=6G
#PBS -o ${job_output_dir}   
#PBS -e ${job_output_dir}

#PBS -m e
#PBS -M ${job_email_recipient}

dataprefix=${final_data_location}
#cd ${dataprefix}/${dataloc}

cd ${seqdir}

NPROCS=16

#if [ "x$PBS_O_WORKDIR" != "x" ]; then
#        cd $PBS_O_WORKDIR
#fi
file=${trim_file}0${PBS_ARRAYID}
cmd="sh ~/atropos_filter_parallel.sh $dir ${file} ${NPROCS}"

date
echo "Executing: $cmd"
eval $cmd
date

