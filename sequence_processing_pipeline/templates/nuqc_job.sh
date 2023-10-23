#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH -p {{queue_name}}
# wall-time-limit in minutes
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}gb
#SBATCH -N {{node_count}}
#SBATCH -c {{cores_per_task}}

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Not operating within an array"
    exit 1
fi

if [[ "${SLURM_ARRAY_TASK_MIN}" -ne 1 ]]; then
    echo "Min array ID is not 1"
    exit 1
fi

if [[ -z ${MMI} ]]; then
    echo "MMI is not set"
    exit 1
fi

if [[ -z ${PREFIX} ]]; then
    echo "PREFIX is not set"
    exit 1
fi

if [[ -z ${OUTPUT} ]]; then
    echo "OUTPUT is not set"
    exit 1
fi

if [[ -z ${TMPDIR} ]]; then
    echo "TMPDIR is not set"
    exit 1
fi

conda activate human-depletion

set -x
set -e

date
hostname
echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}
# output_path = output_path passed to Job objects + 'NuQCJob'
# e.g.: working-directory/ConvertJob, working-directory/QCJob...
cd {{output_path}}

# set a temp directory, make a new unique one under it and
# make sure we clean up as we're dumping to shm
# DO NOT do this casually. Only do a clean up like this if
# you know for sure TMPDIR is what you want.

mkdir -p ${TMPDIR}
export TMPDIR=${TMPDIR}
export TMPDIR=$(mktemp -d)
echo $TMPDIR

function cleanup {
  echo "Removing $TMPDIR"
  rm -fr $TMPDIR
  unset TMPDIR
}
trap cleanup EXIT

export FILES=$(pwd)/$(printf "%s-%d" ${PREFIX} ${SLURM_ARRAY_TASK_ID})
if [[ ! -f ${FILES} ]]; then
    logger ${FILES} not found
    exit 1
fi

delimiter=::MUX::
n=$(wc -l ${FILES} | cut -f 1 -d" ")

for i in $(seq 1 ${n})
do
    line=$(head -n ${i} ${FILES} | tail -n 1)
    r1=$(echo ${line} | cut -f 1 -d" ")
    r2=$(echo ${line} | cut -f 2 -d" ")
    base=$(echo ${line} | cut -f 3 -d" ")
    r1_name=$(basename ${r1} .fastq.gz)
    r2_name=$(basename ${r2} .fastq.gz)

    echo "${i}	${r1_name}	${r2_name}	${base}" >> ${TMPDIR}/id_map

    fastp \
        -l 45 \
        -i ${r1} \
        -I ${r2} \
        -w 7 \
        --adapter_fasta {{knwn_adpt_path}} \
        --html /dev/null \
        --json /dev/null \
        --stdout | \
            sed -r "1~4s/^@(.*)/@${i}${delimiter}\1/"
done > ${TMPDIR}/seqs.fastq

function minimap2 () {
    mmi=$1
    
    echo "$(date) :: $(basename ${mmi})"
    minimap2 -2 -ax sr -t 7 ${mmi} ${TMPDIR}/seqs.fastq | \
        samtools fastq -@ 1 -f 12 -F 256 > ${TMPDIR}/seqs_new.fastq
    mv ${TMPDIR}/seqs_new.fastq ${TMPDIR}/seqs.fastq
}

function runner () {
    i=${1}

    # with the current proposed resources, we likely do not get
    # benefit for parallel gzip write
    mgscripts demux \
        --id-map ${TMPDIR}/id_map \
        --infile ${TMPDIR}/seqs.fastq \
        --output ${OUTPUT} \
        --encoded-id ${i} \
        --threads 1
}
export -f runner

if [[ -f ${MMI} ]]; then
    minimap2 ${MMI}
else
    for mmi in ${MMI}/*.mmi
    do
        minimap2 ${mmi}
    done
fi

mkdir -p ${OUTPUT}

jobs=${SLURM_JOB_CPUS_PER_NODE}

echo "$(date) :: demux start"
# let it do its thing
seq 1 ${n} | parallel -j ${jobs} runner
echo "$(date) :: demux stop"

touch ${OUTPUT}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.completed
