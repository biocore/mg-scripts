#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -p {{queue_name}}
#SBATCH -N {{node_count}}
#SBATCH -n {{nprocs}}
#SBATCH --time {{wall_time_limit}}
#SBATCH --mail-type=ALL
#SBATCH --mail-user qiita.help@gmail.com
#SBATCH --mem-per-cpu {{mem_per_cpu}}
set -x
date
hostname
cd {{run_dir}}
{% if modules_to_load %}
    module load {{modules_to_load}}
{% endif %}
{{cmd_line}}