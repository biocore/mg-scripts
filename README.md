# Sequencing Process Pipeline
Python packaged implementation of our internal pipeline for generating fastqc output from sequence bcl.

run unittests:
cd sequence_processing_pipeline
python -m unittest discover tests

[ToDo: add explanation, installation instructions]

These need to be available to the job-script for human-filtering.
without them, the job won't work. the location of the .mmi file
needs to be hardcoded in fpmmp.beautified.sh. or change the .sh file to
accept a parameter to the mmi file. Also, the .sh file itself
needs to be known to the job object.
scripts/fpmmp.beautified.sh
scripts/human-phix-db.mmi

# Legacy Shell Scripts
Earlier Knight Lab internal Metagenomic processing scripts for demultiplexing, QC and host removal.

The default processing steps are:
1. seq_proc_share.sh: [ToDo: add explanation]
  1.1. bcl_template.sh: [ToDo: add explanation]
  1.2. bowtie_trim_qsub.sh: [ToDo: add explanation]
    1.2.1. filter_job_parallel.sh: [ToDo: add explanation]
  1.3. atropos_filter_parallel.sh: [ToDo: add explanation]
