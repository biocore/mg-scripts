#!/bin/bash -l
#SBATCH -J {{job_name}}             # cleanup
#SBATCH --time {{wall_time_limit}}  # 24:00:00
#SBATCH --mem {{mem_in_gb}}G        # 8G
#SBATCH -N {{node_count}}           # 1
#SBATCH -c {{cores_per_task}}       # 1
#SBATCH -p {{queue_name}}           # qiita

#SBATCH --output tellread-cleanup_%x-%A.out
#SBATCH --error tellread-cleanup_%x-%A.err

if [[ -z "${OUTPUT}" ]]; then
    echo "OUTPUT is not specified"
    exit 1
fi

# remove unused large outputs
rm -rf ${OUTPUT}/biosample_format ${OUTPUT}/1_demult ${OUTPUT}/Full
