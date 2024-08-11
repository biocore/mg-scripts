#!/bin/bash -l

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH --partition=short
#SBATCH -t 96:00:00
#SBATCH -J tellread
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err
#SBATCH --mail-user=qiita.help@gmail.com
#SBATCH --mail-type=BEGIN,FAIL
#SBATCH -p qiita

function logger () { 
    echo "$(date) :: ${@}"; 
    echo "$(date) :: ${@}" 1>&2; 
}

set -x 

if [[ -z "${N_SAMPLES}" ]]; then
    echo "N_SAMPLES is not specified"
    exit 1
fi

if [[ -z "${SEQRUNPATH}" ]]; then
    echo "SEQRUNPATH is not specified"
    exit 1
fi

if [[ -z "${LANE}" ]]; then
    echo "LANE is not specified"
    exit 1
fi

if [[ -z "${SAMPLES}" ]]; then
    echo "SAMPLES is not specified"
    exit 1
fi

if [[ -z "${REFS}" ]]; then
    echo "REFS is not specified"
    exit 1
fi

if [[ -z "${OUTPUT}" ]]; then
    echo "OUTPUT is not specified"
    exit 1
fi

export TMPDIR=/panfs/${USER}/tmp
mkdir -p ${TMPDIR}
export TMPDIR=$(mktemp -d)
seqrun_path=${SEQRUNPATH}

if [[ ${LANE} == "L001" ]]; then
    lane=s_1
elif [[ ${LANE} == "L002" ]]; then
    lane=s_2
elif [[ ${LANE} == "L003" ]]; then
    lane=s_3
elif [[ ${LANE} == "L004" ]]; then
    lane=s_4
elif [[ ${LANE} == "L005" ]]; then
    lane=s_5
elif [[ ${LANE} == "L006" ]]; then
    lane=s_6
elif [[ ${LANE} == "L007" ]]; then
    lane=s_7
elif [[ ${LANE} == "L008" ]]; then
    lane=s_8
else
    echo "Unrecognized lane: ${LANE}"
    exit 1
fi

# yes, hard coded, not great but progress.
extra=""
if [[ ! -z ${REFBASE} ]]; then
    extra="-f ${REFBASE}"
fi

mkdir -p ${OUTPUT}
    
module load singularity_3.6.4
$HOME/qiita-spots/tellread-release-novaseqX/run_tellread_sing.sh \
    -i ${seqrun_path} \
    -o ${OUTPUT} \
    -s $(echo ${SAMPLES} | tr -d '"') \
    -g $(echo ${REFS} | tr -d '"') \
    -j ${SLURM_JOB_CPUS_PER_NODE} \
    ${extra} \
    -l ${lane}

    
if [[ -d ${OUTPUT}/Full ]]; then
    echo "Run appears successful"
elif [[ -d ${OUTPUT}/1_demult/Full ]]; then
    echo "Run appears unsuccessful but has output"
    exit 1
else
    echo "Run appears unsuccessful"
    exit 1
fi
