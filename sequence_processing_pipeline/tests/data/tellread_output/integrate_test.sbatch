#!/bin/bash -l
#SBATCH -J integrate             # integrate
#SBATCH --time 96:00:00  # 24:00:00
#SBATCH --mem 16G        # 8G
#SBATCH -N 1           # 1
#SBATCH -c 4       # 1
#SBATCH -p qiita           # qiita

#SBATCH --output integrate_%x-%A_%a.out
#SBATCH --error integrate_%x-%A_%a.err

# NB SLURM_ARRAY_TASK_ID is exported by Slurm
if [[ -z ${SLURM_ARRAY_TASK_ID} ]]; then
    echo "Not operating in an array"
    exit 1
fi

# NB SLURM_ARRAY_TASK_MIN is exported by Slurm
if [[ ${SLURM_ARRAY_TASK_MIN} -eq 0 ]]; then
    echo "Line extraction assumes 1-based index"
    exit 1
fi

set -x 
set -e
set -o pipefail

samples=($(cat sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/sample_index_list_output.txt | cut -f 2))
sample=${samples[$((${SLURM_ARRAY_TASK_ID} - 1))]} 

# NB TMPDIR IS CREATED IN CURRENT DIRECTORY. CURRENT DIRECTORY MUST BE CORRECT.
export TMPDIR=$(mktemp -d)
function cleanup {                                                              
  echo "Removing $TMPDIR"                                                          
  rm  -r $TMPDIR                                                                   
  unset TMPDIR                                                                  
}                                                                               
trap cleanup EXIT

files=${TMPDIR}/integration.files
/bin/ls -1 sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/Full/*corrected.err_barcode_removed.fastq > ${files}
mkdir -p sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated

if [[ $(grep -c "_R1_${sample}" ${files}) -ne 1 ]]; then
    echo "Multiple matches for ${sample} R1"
    exit 1
fi

if [[ $(grep -c "_R2_${sample}" ${files}) -ne 1 ]]; then
    echo "Multiple matches for ${sample} R2"
    exit 1
fi

if [[ $(grep -c "_I1_${sample}" ${files}) -ne 1 ]]; then
    echo "Multiple matches for ${sample} I1"
    exit 1
fi

r1=$(grep -m 1 "_R1_${sample}" ${files})
r2=$(grep -m 1 "_R2_${sample}" ${files})
i1=$(grep -m 1 "_I1_${sample}" ${files})
r1out=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated/${sample}.R1.fastq.gz
r2out=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated/${sample}.R2.fastq.gz
i1out=sequence_processing_pipeline/tests/2caa8226-cf69-45a3-bd40-1e90ec3d18d0/TRIntegrateJob/integrated/${sample}.I1.fastq.gz

if [[ ! -s ${r1} ]]; then
    echo "${r1} is empty, cannot integrate"
    if [[ -s ${r2} ]]; then
        echo "R1 and R2 are inconsistent"
        exit 1
    fi
    if [[ -s ${i1} ]]; then
        echo "R1 and I1 are inconsistent"
        exit 1
    fi

    # reflect the empties so Qiita can know of them
    touch ${r1out}
    touch ${r2out}
    touch ${i1out}
    exit 0
fi

# this can probably be backgrounded but then you have to get creative to
# not mask a nonzero exit status (e.g., the python process raising)
cat ${i1} | gzip > ${i1out} 

conda activate qp-knight-lab-processing-2022.03
python hello integrate \
    --no-sort \
    --r1-in ${r1} \
    --r2-in ${r2} \
    --i1-in ${i1} \
    --r1-out ${r1out} \
    --r2-out ${r2out} \
    --threads ${SLURM_CPUS_PER_TASK}