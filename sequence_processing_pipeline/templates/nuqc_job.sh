#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH --time {{wall_time_limit_in_min}}
#SBATCH --mem {{mem_in_gb}}gb
#SBATCH -N {{node_count}}
#SBATCH -c {{cores_per_task}}
#SBATCH --output %x-%A_%a.out
#SBATCH --error %x-%A_%a.err

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Not operating within an array"
    exit 1
fi

mamba activate human-depletion

set -x
set -e

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
        --adapter_fasta {{known_adapters_path}} \
        --html /dev/null \
        --json /dev/null \
        --stdout | \
            sed -r "1~4s/^@(.*)/@${i}${delimiter}\1/"
done > ${TMPDIR}/seqs.fastq

for mmi in ${MMI}/*.mmi
do
    echo "$(date) :: $(basename ${mmi})"
    minimap2 -2 -ax sr -t 7 ${mmi} ${TMPDIR}/seqs.fastq | \
        samtools fastq -@ 1 -f 12 -F 256 > ${TMPDIR}/seqs_new.fastq
    echo $(du -sh ${TMPDIR})
    mv ${TMPDIR}/seqs_new.fastq ${TMPDIR}/seqs.fastq
done

mkdir -p ${OUTPUT}

function runner () {
    i=${1}
    python demux-inline.py ${TMPDIR}/id_map ${TMPDIR}/seqs.fastq ${OUTPUT} ${i}
}

export -f runner
jobs=${SLURM_JOB_CPUS_PER_NODE}

echo "$(date) :: demux start"
# let it do its thing
seq 1 ${n} | parallel -j ${jobs} runner
echo "$(date) :: demux stop"

rm -fr ${TMPDIR}
