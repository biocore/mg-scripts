#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH -p {{queue_name}}
### wall-time-limit in minutes
#SBATCH --time {{wall_time_limit}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH -N {{node_count}}
### Note cores_per_task maps to fastp & minimap2 thread counts
### as well as sbatch -c. demux threads remains fixed at 1.
### Note -c set to 4 and thread counts set to 7 during testing.
#SBATCH -c {{cores_per_task}}
#SBATCH --gres=node_jobs:1  # use 1/4 processing slots

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

echo "MMI is ${MMI}"

conda activate human-depletion

set -x
set -e
set -o pipefail

date
hostname
echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}
### output_path = output_path passed to Job objects + 'NuQCJob'
### e.g.: working-directory/ConvertJob, working-directory/QCJob...
cd {{output_path}}
mkdir -p {{html_path}}
mkdir -p {{json_path}}

TMPDIR=/dev/shm
if [[ ! -d ${TMPDIR} ]]; then
	echo "Cannot access ${TMPDIR}"
    exit 1
fi                                                      
export TMPDIR=${TMPDIR}                                                         
export TMPDIR=$(mktemp -d)                                                      
echo $TMPDIR                                                                    
                                                                                
function cleanup {                                                              
  echo "Removing $TMPDIR"                                                          
  rm -fr $TMPDIR                                                                
  unset TMPDIR                                                                  
}                                                                               
trap cleanup EXIT

export FILES=$(printf "%s-%d" ${PREFIX} ${SLURM_ARRAY_TASK_ID})
if [[ ! -f ${FILES} ]]; then
    logger ${FILES} not found
    exit 1
fi

delimiter=::MUX::
export r1_tag=/1
export r2_tag=/2

function mux_runner () {
    n=$(wc -l ${FILES} | cut -f 1 -d" ")

    jobd=${TMPDIR}
    id_map=${jobd}/id_map 
    seqs_r1=${jobd}/seqs.r1.fastq
    seqs_r2=${jobd}/seqs.r2.fastq
    r1_filt=${jobd}/seqs.r1.adapter-removed.fastq
    r2_filt=${jobd}/seqs.r2.adapter-removed.fastq
    
    for i in $(seq 1 ${n})
    do
        line=$(head -n ${i} ${FILES} | tail -n 1)
        r1=$(echo ${line} | cut -f 1 -d" ")
        r2=$(echo ${line} | cut -f 2 -d" ")
        base=$(echo ${line} | cut -f 3 -d" ")
        r1_name=$(basename ${r1} .fastq.gz)
        r2_name=$(basename ${r2} .fastq.gz)

		s_name=$(basename "${r1}" | sed -r 's/\.fastq\.gz//')
		html_name=$(echo "$s_name.html")
		json_name=$(echo "$s_name.json")

        echo "${i}	${r1_name}	${r2_name}	${base}" >> ${id_map}

        fastp \
			-l {{length_limit}} \
			-i ${r1} \
			-I ${r2} \
			-w {{cores_per_task}} \
			--adapter_fasta {{knwn_adpt_path}} \
			--html {{html_path}}/${html_name} \
			--json {{json_path}}/${json_name} \
            --out1 ${r1_filt} \
            --out2 ${r2_filt}

        cat ${r1_filt} | \
            sed -r "1~4s/^@(.*)/@${i}${delimiter}\1/" \
            >> ${seqs_r1} &
        cat ${r2_filt} | \
            sed -r "1~4s/^@(.*)/@${i}${delimiter}\1/" \
            >> ${seqs_r2} &
        wait

        rm ${r1_filt} &
        rm ${r2_filt} &
        wait
    done
}
export -f mux_runner

function minimap2_runner () {
    mmi=$1
    
    echo "$(date) :: $(basename ${mmi})"
    minimap2 -2 -ax sr -t {{cores_per_task}} ${mmi} ${TMPDIR}/seqs.fastq | \
        samtools fastq -@ 1 -f 12 -F 256 > ${TMPDIR}/seqs_new.fastq
    mv ${TMPDIR}/seqs_new.fastq ${TMPDIR}/seqs.fastq

	fcurrent=${1}
    
    jobd=${TMPDIR}
    
    id_map=${jobd}/id_map
    if [[ ! -f ${id_map} ]]; then
        echo "No samples..."
        return
    fi

    seqs_r1=${jobd}/seqs.r1.fastq
    seqs_r2=${jobd}/seqs.r2.fastq
    seqs_mmpe_r1=${jobd}/seqs.mmpe.r1.fastq
    seqs_mmpe_r2=${jobd}/seqs.mmpe.r2.fastq
    seqs_mmpese=${jobd}/seqs.mmpese.fastq
    seqs_mmpese_r1=${jobd}/seqs.mmpese.r1.fastq
    seqs_mmpese_r2=${jobd}/seqs.mmpese.r2.fastq

    # PE operation
    minimap2 -2 -ax sr -t 12 \
        ${fcurrent} \
        ${seqs_r1} \
        ${seqs_r2} | \
            samtools fastq \
                -@ 4 \
                -f 12 \
                -F 256 \
                -N \
                -1 ${seqs_mmpe_r1} \
                -2 ${seqs_mmpe_r2}

    rm ${seqs_r1} &
    rm ${seqs_r2} &
    # no need to block

    # SE operation
    # run r1/r2 serially to avoid minimap2 picking up on interleve
    minimap2 -2 -ax sr -t 12 \
        ${fcurrent} \
        <(cat ${seqs_mmpe_r1} ${seqs_mmpe_r2}) | \
            samtools fastq \
                -@ 4 \
                -f 4 \
                -F 256 \
                -0 ${seqs_mmpese}

    rm ${seqs_mmpe_r1} &
    rm ${seqs_mmpe_r2} &
    # no need to block

    ./splitter/target/release/splitter ${seqs_mmpese} ${seqs_mmpese_r1} ${delimiter} ${r1_tag} &
    ./splitter/target/release/splitter ${seqs_mmpese} ${seqs_mmpese_r2} ${delimiter} ${r2_tag} &
    wait

    rm ${seqs_mmpese} &

    fastq_pair -t 50000000 ${seqs_mmpese_r1} ${seqs_mmpese_r2}
    rm ${seqs_mmpese_r1}.single.fq &
    rm ${seqs_mmpese_r2}.single.fq &
    rm ${seqs_mmpese_r1} &
    rm ${seqs_mmpese_r2} &
    wait

    mv ${seqs_mmpese_r1}.paired.fq ${seqs_r1}
    mv ${seqs_mmpese_r2}.paired.fq ${seqs_r2}
}
export -f minimap2_runner

function demux_runner () {
    n_demux_jobs=${SLURM_CPUS_PER_TASK}
    jobd=${TMPDIR}
    id_map=${jobd}/id_map
    seqs_r1=${jobd}/seqs.r1.fastq
    seqs_r2=${jobd}/seqs.r2.fastq

    id_map=${jobd}/id_map
    if [[ ! -f ${id_map} ]]; then
        echo "No samples..."
        return
    fi

    for idx in $(seq 0 ${n_demux_jobs})
    do
        python {{demux_path}} \
            ${id_map} \
            <(cat ${seqs_r1} ${seqs_r2}) \
            ${OUTPUT} \
            ${idx} \
            ${n_demux_jobs} &
    done
    wait
}
export -f demux_runner

mmi_files=${TMPDIR}/mmi-files
if [[ -d ${MMI} ]]; then
    /bin/ls -1 ${MMI}/*.mmi > ${mmi_files}
else
    echo ${MMI} > ${mmi_files}
fi
n_files=$(wc -l ${mmi_files} | cut -f 1 -d" ")

mux-runner

# process with minimap2
for idx in $(seq 1 ${n_files})
do
    fcurrent=$(head -n ${idx} ${mmi_files} | tail -n 1)
    export fcurrent
    echo "$(date) :: $(basename ${fcurrent})"
    minimap2-runner
done

mkdir -p ${OUTPUT}

echo "$(date) :: demux start"
demux-runner
echo "$(date) :: demux stop"

touch ${OUTPUT}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.completed
