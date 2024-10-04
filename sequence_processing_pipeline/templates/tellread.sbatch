#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH -p {{queue_name}}
#SBATCH -N {{node_count}}
#SBATCH -c {{cores_per_task}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH --time {{wall_time_limit}}

#SBATCH --output tellread_%x-%A.out
#SBATCH --error tellread_%x-%A.err

set -x

export TMPDIR={{tmp_dir}}
mkdir -p ${TMPDIR}
export TMPDIR=$(mktemp -d)

mkdir -p {{output}}

module load {{modules_to_load}}
{{sing_script_path}} \
    -i {{rundir_path}} \
    -o {{output}} \
    -s $(echo {{samples}} | tr -d '"') \
    -g $(echo {{refs}} | tr -d '"') \
    -j ${SLURM_JOB_CPUS_PER_NODE} {{extra}} \
    -l {{lane}}

if [[ -d {{output}}/Full ]]; then
    echo "Run appears successful"
elif [[ -d {{output}}/1_demult/Full ]]; then
    echo "Run appears unsuccessful but has output"
    exit 1
else
    echo "Run appears unsuccessful"
    exit 1
fi
