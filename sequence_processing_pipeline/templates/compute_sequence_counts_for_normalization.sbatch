#!/bin/bash -l
#SBATCH -J {{job_name}}             # norm
#SBATCH --time {{wall_time_limit}}  # 24:00:00
#SBATCH --mem {{mem_in_gb}}G        # 8G
#SBATCH -N {{node_count}}           # 1
#SBATCH -c {{cores_per_task}}       # 1
#SBATCH -p {{queue_name}}           # qiita

#SBATCH --output {{output}}/logs/compute_sequence_counts_%x-%A_%a.out
#SBATCH --error {{output}}/logs/compute_sequence_counts_%x-%A_%a.err

# NB: output appears normal w/out.
# source activate qiime2-2023.5

set -x
set -e
set -o pipefail

echo $TMPDIR

mkdir -p {{output_path}}
wc -l {{tellread_output}}/Full/*_I1_C5[0-9][0-9].fastq.gz.corrected.err_barcode_removed.fastq > {{output_path}}/record_counts.txt
python {{plot_counts_path}} {{output_path}}/record_counts.txt {{sample_sheet}} {{output_path}}

conda activate qp-knight-lab-processing-2022.03
python {{create_picklist_path}} {{read_counts_path}}
