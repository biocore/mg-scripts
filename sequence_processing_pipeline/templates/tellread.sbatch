#!/bin/bash -l
#SBATCH -J {{job_name}}
#SBATCH -p {{queue_name}}
#SBATCH -N {{node_count}}
#SBATCH -c {{cores_per_task}}
#SBATCH --mem {{mem_in_gb}}G
#SBATCH --time {{wall_time_limit}}

#SBATCH --output {{output}}/logs/tellread_%x-%A.out
#SBATCH --error {{output}}/logs/tellread_%x-%A.err

set -x

module load {{modules_to_load}}
{{sing_script_path}} \
    -i {{rundir_path}} \
    -o {{output}} \
    -s $(echo {{samples}} | tr -d '"') \
    -g $(echo {{refs}} | tr -d '"') \
    -j ${SLURM_JOB_CPUS_PER_NODE} {{extra}} \
    -l {{lane}}

if [[ -d {{output}}/Full ]]; then
    echo "tellread.sbatch successful"
else
    echo "tellread.sbatch unsuccessful"
    exit 1
fi
