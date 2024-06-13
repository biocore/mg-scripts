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
#SBATCH --gres=node_jobs:{{gres_value}}


echo "---------------"
echo "Run details:"
echo "$SLURM_JOB_NAME $SLURM_JOB_ID $SLURMD_NODENAME $SLURM_ARRAY_TASK_ID"
echo "---------------"

if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Not operating within an array"
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

conda activate qp-knight-lab-processing-2022.03
module load {{modules_to_load}}

set -x
set -e
set -o pipefail

export FILES=$(printf "%s-%d" ${PREFIX} ${SLURM_ARRAY_TASK_ID})
if [[ ! -f ${FILES} ]]; then
    logger ${FILES} not found
    exit 1
fi
# set a temp directory, make a new unique one under it and
# make sure we clean up as we're dumping to shm
# DO NOT do this casually. Only do a clean up like this if
# you know for sure TMPDIR is what you want.

WKDIR=${OUTPUT}/
TMPDIR=${OUTPUT}
export TMPDIR=${TMPDIR}
export TMPDIR=$(mktemp -d)
echo $TMPDIR

mkdir -p ${WKDIR}NuQCJob/fastp_reports_dir/html
mkdir -p ${WKDIR}NuQCJob/fastp_reports_dir/json

export ADAPTER_ONLY_OUTPUT=${OUTPUT}/only-adapter-filtered
mkdir -p ${ADAPTER_ONLY_OUTPUT}

function cleanup {
  echo "Removing $TMPDIR"
  rm -fr $TMPDIR
  unset TMPDIR
}
trap cleanup EXIT

export delimiter=::MUX::
export r1_tag=/1
export r2_tag=/2
function mux-runner () {
    n=$(wc -l ${FILES} | cut -f 1 -d" ")

    jobd=${TMPDIR}
    id_map=${jobd}/id_map
    seqs_r1=${jobd}/seqs.r1.fastq.gz
    seqs_r2=${jobd}/seqs.r2.fastq
    r1_filt=${jobd}/seqs.r1.adapter-removed.fastq.gz

    for i in $(seq 1 ${n})
    do
        line=$(head -n ${i} ${FILES} | tail -n 1)
        r1=$(echo ${line} | cut -f 1 -d" ")
        r2=$(echo ${line} | cut -f 2 -d" ")
        base=$(echo ${line} | cut -f 3 -d" ")
        r1_name=$(basename ${r1} .fastq.gz)
        r2_name=$(basename ${r2} .fastq.gz)
        r1_adapter_only=${ADAPTER_ONLY_OUTPUT}/${r1_name}.fastq.gz

        s_name=$(basename "${r1}" | sed -r 's/\.fastq\.gz//')
        html_name=$(echo "$s_name.html")
        json_name=$(echo "$s_name.json")

        echo -e "${i}\t${r1_name}\t${r2_name}\t${base}" >> ${id_map}

        fastp \
            -l {{length_limit}} \
            -i ${r1} \
            -I ${r2} \
            -w {{cores_per_task}} \
            --adapter_fasta {{knwn_adpt_path}} \
            --html {{html_path}}/${html_name} \
            --json {{json_path}}/${json_name} \
            --stdout | gzip > ${r1_filt}

        # multiplex and write adapter filtered data all at once
        zcat ${r1_filt} | \
            sed -r "1~4s/^@(.*)/@${i}${delimiter}\1/" \
            >> ${seqs_r1} 
        cat ${r1_filt} | \
            gzip -c > ${r1_adapter_only} &
        wait

        rm ${r1_filt} &
        wait
    done

    # minimap/samtools pair commands are now generated in NuQCJob._generate_mmi_filter_cmds()
    # and passed to this template. This method assumes ${jobd} is the correct location to
    # filter files, the initial file is "${jobd}/seqs.r1.fastq"), and the output name is
    # "${jobd}/seqs.r1.ALIGN.fastq".
    {{mmi_filter_cmds}}

    {{movi_path}} query \
        --index /scratch/movi_hg38_chm13_hprc94 \
        --read ${jobd}/seqs.r1.ALIGN.fastq \
        --stdout > ${jobd}/seqs.movi.txt
        
    python {{pmls_path}} <(zcat ${jobd}/seqs.movi.txt.gz) | \
        seqtk subseq ${jobd}/seqs.r1.ALIGN.fastq.gz - | gzip > ${jobd}/seqs.r1.final.fastq.gz
         
    {{splitter_binary}} ${jobd}/seqs.r1.final.fastq \
        ${jobd}/reads.r1.fastq ${delimiter} ${r1_tag} &
    {{splitter_binary}} ${jobd}/seqs.r1.final.fastq \
        ${jobd}/reads.r2.fastq ${delimiter} ${r2_tag} &
    wait
    fastq_pair -t 50000000 ${jobd}/reads.r1.fastq ${jobd}/reads.r2.fastq

    # keep seqs.movi.txt and migrate it to NuQCJob directory.
    mv ${jobd}/seqs.movi.txt {{output_path}}/seqs.movi.${SLURM_ARRAY_TASK_ID}.txt
}
export -f mux-runner


function demux-runner () {
    n_demux_jobs=${SLURM_CPUS_PER_TASK}
    jobd=${TMPDIR}
    id_map=${jobd}/id_map
    seqs_r1=${jobd}/reads.r1.fastq.paired.fq
    seqs_r2=${jobd}/reads.r2.fastq.paired.fq

    id_map=${jobd}/id_map
    if [[ ! -f ${id_map} ]]; then
        echo "No samples..."
        return
    fi

    for idx in $(seq 0 ${n_demux_jobs})
    do
        python {{demux_path}} \
            --id-map ${id_map} \
            --infile <(cat ${seqs_r1} ${seqs_r2}) \
            --output ${OUTPUT} \
            --task ${idx} \
            --maxtask ${n_demux_jobs} &
    done
    wait
}
export -f demux-runner

mux-runner

mkdir -p ${OUTPUT}

echo "$(date) :: demux start"
demux-runner
echo "$(date) :: demux stop"

touch ${OUTPUT}/${SLURM_JOB_NAME}.${SLURM_ARRAY_TASK_ID}.completed
