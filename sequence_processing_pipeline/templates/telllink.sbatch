#!/bin/bash -l
#SBATCH --mem 160G
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 96:00:00
#SBATCH -J tellink
#SBATCH --output %x-%A_%a.out
#SBATCH --error %x-%A_%a.err
#SBATCH --mail-user=qiita.help@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p qiita

set -x 
set -e

module load singularity_3.6.4

if [[ -z "${LABELTAG}" ]]; then
    echo "LABEL is not specified"
    exit 1
fi

base=/panfs/${USER}/${LABELTAG}
if [[ ! -d ${base} ]]; then
    echo "${base} not found"
    exit 1
fi

samples=($(cat ${base}/sample_index_list_${LABELTAG}.txt | cut -f 2))
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]} 

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

HOME_PATH=/projects/long_read_collab/code/tellseq/release_v1.11/
${HOME_PATH}/tellink-release/run_tellink_sing.sh \
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

