# Sequencing Process Pipeline
Python packaged implementation of our internal pipeline for generating fastqc output from sequence bcl.

[ToDo: add explanation, installation instructions]

# Legacy Shell Scripts
Earlier Knight Lab internal Metagenomic processing scripts for demultiplexing, QC and host removal.

The default processing steps are:
1. seq_proc_share.sh: [ToDo: add explanation]
  1.1. bcl_template.sh: [ToDo: add explanation]
  1.2. bowtie_trim_qsub.sh: [ToDo: add explanation]
    1.2.1. filter_job_parallel.sh: [ToDo: add explanation]
  1.3. atropos_filter_parallel.sh: [ToDo: add explanation]
