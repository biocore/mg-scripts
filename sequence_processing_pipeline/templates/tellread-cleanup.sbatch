#!/bin/bash -l
#SBATCH -J {{job_name}}             # cleanup
#SBATCH --time {{wall_time_limit}}  # 24:00:00
#SBATCH --mem {{mem_in_gb}}G        # 8G
#SBATCH -N {{node_count}}           # 1
#SBATCH -c {{cores_per_task}}       # 1
#SBATCH -p {{queue_name}}           # qiita

# for now comment these out as qiita is responsible for notifying users.
###SBATCH --mail-user=qiita.help@gmail.com
###SBATCH --mail-type=BEGIN,FAIL

# for now these can be left hard-coded.
#SBATCH --output %x-%A.out
#SBATCH --error %x-%A.err

if [[ -z "${OUTPUT}" ]]; then
    echo "OUTPUT is not specified"
    exit 1
fi

# remove unused large outputs
rm -rf ${OUTPUT}/biosample_format ${OUTPUT}/1_demult ${OUTPUT}/Full
