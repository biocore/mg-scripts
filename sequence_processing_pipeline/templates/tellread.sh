#!/bin/bash

script_name=${0##*/}

function help () {
    echo "Submit for TELL-read"
    echo ""
    echo "Usage: ${script_name} -s <seqrunpath> -l <lane> [-r reference_map] [-b reference_base]"
    echo ""
    echo -e "\t-s\tPath to the sequencing run."  
    echo -e "\t-i\tThe sample sheet."  
    echo -e "\t-l\tThe lane to process."  
    echo -e "\t-r\tA file specifying reference genomes to use [OPTIONAL]"  
    echo -e "\t-b\tReference genome base directory [OPTIONAL]"  
    echo -e "\t-m\tMode, isolate or metagenomic [OPTIONAL]"  
    echo "" 
}

# references right now are only used for techdev

# derived from https://www.redhat.com/sysadmin/arguments-options-bash-scripts
while getopts "hs:i:l:r:b:m:" option; do
   case ${option} in
    h) 
        help
        exit;;
    s) seqrunpath=${OPTARG};;
    l) lane=${OPTARG};;
    r) reference_map=${OPTARG};;
    b) reference_base=${OPTARG};;
    m) mode=${OPTARG};;
    \?)
         echo "Error: Invalid option"
         exit;;
    *)
         echo "Error: Invalid option"
         exit;;
   esac
done

# nifty
# https://unix.stackexchange.com/a/621007
: ${seqrunpath:?Missing -s}
: ${lane:?Missing -i}

if [[ ! -z ${reference_map} || ! -z ${reference_base} ]]; then
    if [[ -z ${reference_map} ]]; then
        echo "-b used without -r"
        exit 1
    fi
    if [[ -z ${reference_base} ]]; then
        echo "-r used without -b"
        exit 1
    fi
    if [[ ! -d ${reference_base} ]]; then
        echo "reference base not found"
        exit 1
    fi

    tag=reference-based
else
    tag=reference-free
fi

samplesheet="/home/qiita_test/qiita-spots/tellread_mapping.csv"

# trim trailing slash
# https://stackoverflow.com/a/32845647/19741
safepath=$(echo ${seqrunpath} | sed 's:/*$::')  
label=$(basename ${safepath})
labeltag=${label}-${tag}
output=/panfs/${USER}/${labeltag}

if [[ ! -d ${seqrunpath}/Data/Intensities/BaseCalls/${lane} ]]; then
    echo "Cannot access the lane"
    exit 1
fi

if [[ ${seqrunpath} == *"_iSeq_Runs"* ]]; then
    echo "FOO"
    sbatch_cores=2
    sbatch_mem=8G
    norm=TRUE
    wall=24:00:00
    mode=NA
elif [[ ${seqrunpath} == *"_MiSeq_Runs"* ]]; then
    echo "BAR"
    sbatch_cores=2
    sbatch_mem=8G
    norm=TRUE
    wall=24:00:00
    mode=NA
else
    echo "BAZ"
    sbatch_cores=16
    sbatch_mem=160G
    norm=FALSE
    assemble=TRUE
    wall=48:00:00
fi

if [[ ${mode} == "isolate" ]]; then
    ISOLATE_MODE=TRUE
elif [[ ${mode} == "metagenomic" ]]; then
    ISOLATE_MODE=FALSE
elif [[ ${mode} == "NA" ]]; then
    ISOLATE_MODE=FALSE
else
    echo "unknown mode: ${mode}"
    exit 1
fi

set -e
set -o pipefail

declare -a s
declare -a g
# below extended regex might be broken because C5\d\d happens in column 0, not column 1
# of the hacked sample-sheet. 
#for sample in $(egrep -o ",C5[0-9][0-9]," ${samplesheet} | tr -d "," | sort)
for sample in $(egrep -o "^C5.*," ${samplesheet} | tr -d "," | sort)
do
    echo "SAMPLE: ${sample}"
    # get references if they exist
    if [[ -f ${reference_map} ]]; then
        if $(grep -Fq ${sample} ${reference_map}); then
            ref=$(grep -m 1 ${sample} ${reference_map} | cut -f 2 -d"," | tr -d "\n")
            if [[ ${ref} != "NONE" ]]; then
                if [[ ! -d "${reference_base}/${ref}" ]]; then
                    echo "${reference_base}/${ref}"
                    echo "${ref} not found"
                    exit 1
                fi
                g[${#g[@]}]=${ref}
                s[${#s[@]}]=${sample}
            fi
        fi
    else
        g[${#g[@]}]=NONE
        s[${#s[@]}]=${sample}
    fi
done
n_samples=${#s[@]}

echo "Submitting:"
echo "S: ${s[@]}"
echo "G: ${g[@]}"

# https://stackoverflow.com/a/17841619/19741
function join_by { local IFS="$1"; shift; echo "$*"; }
s=$(join_by , "${s[@]}")
g=$(join_by , "${g[@]}")

base=$(dirname ${0})
submit_script=$(dirname ${0})/tellread.sbatch
integrate_script=$(dirname ${0})/integrate.sbatch
norm_script=$(dirname ${0})/compute_sequence_counts_for_normalization.sbatch
asm_cloudspades_script=$(dirname ${0})/cloudspades.sbatch
clean_script=$(dirname ${0})/tellread-cleanup.sbatch

if [[ ${ISOLATE_MODE} == "TRUE" ]]; then
    asm_tellink_script=$(dirname ${0})/telllink-isolate.sbatch
    asm_cloudspades_script=$(dirname ${0})/cloudspades-isolate.sbatch
else
    asm_cloudspades_script=$(dirname ${0})/cloudspades.sbatch
    asm_tellink_script=$(dirname ${0})/telllink.sbatch
fi

if [[ ! -f ${submit_script} ]]; then
    echo "Cannot access submit script"
    exit 1
fi
if [[ ! -f ${asm_cloudspades_script} ]]; then
    echo "Cannot access cloudspades assembly script"
    exit 1
fi
if [[ ! -f ${asm_tellink_script} ]]; then
    echo "Cannot access tell-link assembly script"
    exit 1
fi
if [[ ! -f ${integrate_script} ]]; then
    echo "Cannot access integrate script"
    exit 1
fi
if [[ ! -f ${clean_script} ]]; then
    echo "Cannot access clean script"
    exit 1
fi

datetag=$(date "+%Y.%m.%d")
scriptcopy=$(pwd)/tellread_script-${datetag}.sh
submitcopy=$(pwd)/tellread_submission-${datetag}.sbatch
asmcscopy=$(pwd)/assembly_submission_cloudspades-${datetag}.sbatch
asmtlcopy=$(pwd)/assembly_submission_tell-link-${datetag}.sbatch
normcopy=$(pwd)/norm_submission-${datetag}.sbatch
intcopy=$(pwd)/integrate_submission-${datetag}.sbatch
cleancopy=$(pwd)/tellread-cleanup-${datetag}.sbatch
arguments=$(pwd)/provided_script_arguments.txt
if [[ -f ${scriptcopy} ]]; then
    echo "Existing script copy ${scriptcopy} found, not overwriting, delete to resubmit"
    exit 1
fi
if [[ -f ${submitcopy} ]]; then
    echo "Existing submission ${submitcopy} found, not overwriting, delete to resubmit"
    exit 1
fi

# CHARLIE
echo $@ > ${arguments}
cp ${0} ${scriptcopy}
cp ${submit_script} ${submitcopy}
cp ${asm_cloudspades_script} ${asmcscopy}
cp ${asm_tellink_script} ${asmtlcopy}
cp ${integrate_script} ${intcopy}
cp ${clean_script} ${cleancopy}
chmod gou-w ${scriptcopy} ${submitcopy} ${asmcopy} ${intcopy} ${arguments} ${cleancopy}

set -x

echo "C"

trjob=$(sbatch \
          --parsable \
          -J ${labeltag}-${datetag} \
          -c ${sbatch_cores} \
          --mem ${sbatch_mem} \
          --time ${wall} \
          --export BASE=${base},N_SAMPLES=${n_samples},SEQRUNPATH=${seqrunpath},LANE=${lane},REFMAP=${reference_map},REFBASE=${reference_base},OUTPUT=${output},SAMPLES=\"${s}\",REFS=\"${g}\" \
          ${submit_script})

echo "D"

if [[ ${norm} == "TRUE" ]]; then
    cp ${norm_script} ${normcopy}
    chmod gou-w ${normcopy}
    norm_counts_job=$(sbatch \
                        --parsable \
                        --dependency=afterok:${trjob} \
                        -J ${labeltag}-${datetag}-norm-counts \
                        --export BASE=${base},TELLREAD_OUTPUT=${output},OUTPUT=$(pwd),SAMPLESHEET=${samplesheet} \
                        ${norm_script})
fi

echo "E"
integrate_job=$(sbatch \
                    --parsable \
                    -J ${labeltag}-${datetag}-integrate \
                    --dependency=afterok:${trjob} \
                    --array 1-${n_samples} \
                    --export BASE=${base},LABELTAG=${labeltag},OUTPUT=${output} \
                    ${integrate_script})

if [[ ${assemble} == "TRUE" ]]; then
    csj=$(sbatch \
            --parsable \
            --dependency=aftercorr:${integrate_job} \
            -J ${labeltag}-${datetag}-cloudspades \
            --array 1-${n_samples} \
            --export LABELTAG=${labeltag},OUTPUT=${output} \
            ${asm_cloudspades_script})
    tlj=$(sbatch \
            --parsable \
            --dependency=aftercorr:${integrate_job} \
            -J ${labeltag}-${datetag}-tell-link \
            --array 1-${n_samples} \
            --export LABELTAG=${labeltag},OUTPUT=${output} \
            ${asm_tellink_script})
    cleanupdep=${csj}:${tlj}
else
    cleanupdep=${integrate_job}
    echo "Not assembling"
fi

cleanup=$(sbatch \
            --parsable \
            -J ${labeltag}-${datetag}-cleanup \
            --dependency=afterok:${cleanupdep} \
            --export OUTPUT=${output} \
            ${clean_script})
