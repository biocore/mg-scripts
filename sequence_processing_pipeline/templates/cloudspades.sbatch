#!/bin/bash -l
#SBATCH -J {{job_name}}             # cs-assemble
#SBATCH --time {{wall_time_limit}}  # 24:00:00
#SBATCH --mem {{mem_in_gb}}G        # 128G
#SBATCH -N {{node_count}}           # 1
#SBATCH -c {{cores_per_task}}       # 12
#SBATCH -p {{queue_name}}           # qiita

#SBATCH --output cloudspades_%x-%A_%a.out
#SBATCH --error cloudspades_%x-%A_%a.err

source activate qiime2-2023.5

set -x 
set -e

module load {{modules_to_load}}

samples=($(cat {{output_path}}/sample_index_list_output.txt | cut -f 2))            

# assumes 1-based array index, eg --array 1-N
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]} 

cs={{output_path}}/cloudspades/${sample}

if [[ ! -z ${FORCE} && ${FORCE} == "TRUE" ]]; then
    if [[ -d ${cs} ]]; then
        rm -fr ${cs}
    fi
fi

mkdir -p ${cs}
pushd {{cloudspades_path}}/assembler/bin

./spades.py \
    -o ${cs} \
    --gemcode1-1 {{output_path}}/integrated/${sample}.R1.fastq.gz \
    --gemcode1-2 {{output_path}}/integrated/${sample}.R2.fastq.gz \
    --meta \
    -t ${SLURM_JOB_CPUS_PER_NODE} > ${cs}/stdoutstderr.log 2>&1
module unload gcc_9.3.0
popd

# TODO: Look for alternative method to load quast
#mamba activate quast
#quast \
#    -o ${cs}/quast-scaffolds \
#    -t ${SLURM_JOB_CPUS_PER_NODE} \
#    ${cs}/scaffolds.fasta > ${cs}/quast-stdoutstderr.log 2>&1

# remove intermediates that currently dont have a downstream use
#if [[ -d ${cs}/K21 ]]; then
#    rm -fr ${cs}/K21 ${cs}/K33 ${cs}/K55 ${cs}/corrected ${cs}/tmp
#fi
