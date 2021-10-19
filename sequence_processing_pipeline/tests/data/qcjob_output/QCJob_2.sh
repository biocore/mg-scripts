#!/bin/bash
#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_Feist_11661
#PBS -q queue_name
#PBS -l nodes=1:ppn=16
#PBS -V
#PBS -l walltime=24:00:00
#PBS -l mem=8gb
#PBS -o localhost:sequence_processing_pipeline/tests/data/MyRunDir/QCJob_2.out.log.${PBS_ARRAYID}
#PBS -e localhost:sequence_processing_pipeline/tests/data/MyRunDir/QCJob_2.err.log.${PBS_ARRAYID}
#PBS -t 1-108%30
set -x
date
hostname
echo ${PBS_JOBID} ${PBS_ARRAYID}
cd sequence_processing_pipeline/tests/data/MyRunDir
module load 
offset=${PBS_ARRAYID}
step=$(( $offset - 0 ))
cmd0=$(head -n $step sequence_processing_pipeline/tests/data/MyRunDir/split_file_Feist_11661.array-details | tail -n 1)
eval $cmd0