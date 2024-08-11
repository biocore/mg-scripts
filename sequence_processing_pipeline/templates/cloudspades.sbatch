#!/bin/bash -l
#SBATCH -J cs-assemble
#SBATCH --time 24:00:00
#SBATCH --mem 128gb
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --output %x-%A_%a.out
#SBATCH --error %x-%A_%a.err
#SBATCH --mail-user=qiita.help@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p qiita

source activate qiime2-2023.5
function logger () { 
    echo "$(date) :: ${@}"; 
    echo "$(date) :: ${@}" 1>&2; 
}

set -x 
set -e

echo $TMPDIR

if [[ -z "${LABELTAG}" ]]; then
    echo "LABELTAG is not specified"
    exit 1
fi

base=${OUTPUT}
if [[ ! -d ${base} ]]; then
    echo "${base} not found"
    exit 1
fi

mamba activate activate qiime2-2023.5
module load gcc_9.3.0 

samples=($(cat ${base}/sample_index_list_${LABELTAG}.txt | cut -f 2))            

# assumes 1-based array index, eg --array 1-N
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]} 

cs=${base}/cloudspades/${sample}

if [[ ! -z ${FORCE} && ${FORCE} == "TRUE" ]]; then
    if [[ -d ${cs} ]]; then
        rm -fr ${cs}
    fi
fi

mkdir -p ${cs}
pushd /home/qiita_test/qiita-spots/spades-cloudspades-0.1/assembler/bin

./spades.py \
    -o ${cs} \
    --gemcode1-1 ${base}/integrated/${sample}.R1.fastq.gz \
    --gemcode1-2 ${base}/integrated/${sample}.R2.fastq.gz \
    --meta \
    -t ${SLURM_JOB_CPUS_PER_NODE} > ${cs}/stdoutstderr.log 2>&1
module unload gcc_9.3.0
popd

mamba activate quast                                                           
quast \
    -o ${cs}/quast-scaffolds \
    -t ${SLURM_JOB_CPUS_PER_NODE} \
    ${cs}/scaffolds.fasta > ${cs}/quast-stdoutstderr.log 2>&1

# remove intermediates that currently dont have a downstream use
if [[ -d ${cs}/K21 ]]; then
    rm -fr ${cs}/K21 ${cs}/K33 ${cs}/K55 ${cs}/corrected ${cs}/tmp
fi
