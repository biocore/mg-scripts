#!/bin/bash -l
#SBATCH -J {{job_name}}             # tellink
#SBATCH --mem {{mem_in_gb}}G        # 160G
#SBATCH -N {{node_count}}           # 1
#SBATCH -c {{cores_per_task}}       # 16
#SBATCH --time {{wall_time_limit}}  # 96:00:00
#SBATCH -p {{queue_name}}           # qiita

#SBATCH --output telllink_%x-%A_%a.out
#SBATCH --error telllink_%x-%A_%a.err

set -x 
set -e

module load {{modules_to_load}}

if [[ -z "${LABELTAG}" ]]; then
    echo "LABEL is not specified"
    exit 1
fi

base={{output_path}}
if [[ ! -d ${base} ]]; then
    echo "${base} not found"
    exit 1
fi

samples=($(cat ${base}/sample_index_list_output.txt | cut -f 2))
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]} 

# TODO: leave these hardcoded for now
k=79
lc=35
cores=${SLURM_CPUS_PER_TASK}

tl=${base}/tell-link/${sample}
if [[ ! -z ${FORCE} && ${FORCE} == "TRUE" ]]; then
    if [[ -d ${tl} ]]; then
        rm -fr ${tl}
    fi
fi

mkdir -p ${tl}

{{sing_path}} \
    -r1 ${base}/integrated/${sample}.R1.fastq.gz \
    -r2 ${base}/integrated/${sample}.R2.fastq.gz \
    -i1 ${base}/integrated/${sample}.I1.fastq.gz \
    -d metagenomics \
    -o ${tl}/${LABELTAG}-link_${sample}_global_${k}_local_${lc} \
    -k ${k} \
    -lc ${lc} \
    -p ${sample} \
    -j ${cores}

# remove temporary data
if [[ -d ${tl}/${LABELTAG}-link_${sample}_global_${k}_local_${lc}/${sample}/__skipping ]]; then
    rm -fr ${tl}/${LABELTAG}-link_${sample}_global_${k}_local_${lc}/${sample}/__skipping
fi
