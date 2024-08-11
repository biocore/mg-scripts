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

### --gres=gpu:1
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

#acs=${base}/cloudspades/${sample}-ariadne
#acscs=${acs}/assembled
mkdir -p ${cs}
#mkdir -p ${acs}
#mkdir -p ${acscs}

#pushd ~/spades-cloudspades-paper/assembler/
#pushd /home/mcdonadt/cloudspades-0.1/spades-cloudspades-0.1/assembler/bin
#pushd /home/qiita/CHARLIE/TELLREAD/spades-cloudspades-0.1/assembler/bin
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

#pushd $HOME/2023.08.29-ariadne/ariadne
#mamba activate ariadne-gcc8.5.0
#OLD_LDD_LIBRARY_PATH=${LDD_LIBRARY_PATH}
#export LDD_LIBRARY_PATH=${HOME}/miniconda3/envs/ariadne-gcc8.5.0/include/:${LDD_LIBRARY_PATH}
## parameters from Lauren Mak, 9.25.23
#./spades.py \
#    -o ${acs} \
#    --gemcode1-1 ${base}/integrated/${sample}.R1.fastq.gz \
#    --gemcode1-2 ${base}/integrated/${sample}.R2.fastq.gz \
#    --assembly-graph ${cs}/assembly_graph_with_scaffolds.gfa \
#    --meta \
#    --only-assembler \
#    --search-distance 5000 --size-cutoff 6 -k 55 \
#    -t ${SLURM_JOB_CPUS_PER_NODE} > ${acs}/stdoutstderr.log 2>&1
#export LDD_LIBRARY_PATH=${OLD_LDD_LIBRARY_PATH}
#popd
#
#pushd ~/spades-cloudspades-paper/assembler/
#module load gcc_9.3.0 
#./spades.py \
#    -o ${acscs} \
#    --gemcode1-1 ${acs}/K55/5000.R1.fastq \
#    --gemcode1-2 ${acs}/K55/5000.R2.fastq \
#    --meta \
#    -t ${SLURM_JOB_CPUS_PER_NODE} > ${acscs}/stdoutstderr.log 2>&1
#module unload gcc_9.3.0
#popd
#
#mamba activate quast                                                           
#quast \
#    -o ${acscs}/quast-scaffolds \
#    -t ${SLURM_JOB_CPUS_PER_NODE} \
#    ${acscs}/scaffolds.fasta > ${acscs}/quast-stdoutstderr.log 2>&1
