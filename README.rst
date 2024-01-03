# Sequence Processing Pipeline

A Jupyter notebook to assist wet lab shotgun pipeline.
A packaged Python-based implementation of Knight Lab's sequencing process pipeline.

## Installation

To install this package, first clone this repository from GitHub:

```bash
git clone https://github.com/biocore/mg-scripts.git
```

Create a Python3 Conda environment in which to run the notebook:

```bash
conda create -n sp_pipeline 'python==3.9' numpy pandas click scipy matplotlib fastq-pair
```

Activate the Conda environment:

```bash
source activate sp_pipeline
```

Change directory to the cloned repository folder and install:

```bash
cd mg-scripts
pip install -e .
```

This will automatically install https://github.com/biocore/metagenomics_pooling_notebook.git, a dependency of mg-scripts and the sequence_processing_pipeline.

## Running Unittests

Change directory to the downloaded repository folder:

```bash
cd mg-scripts
nosetests --with-coverage --cover-inclusive --cover-package sequence_processing_pipeline
```

## Getting Started

Review Pipeline.py and main.py to learn how to import and access package functionality:

```bash
cd mg-scripts/sequence_processing_pipeline
more Pipeline.py
more main.py
```

Adjust configuration settings as needed:

```bash
cd mg-scripts/sequence_processing_pipeline
vi configuration.json
```

Please note that the setting 'minimap2_databases' is expected to be a list of paths to individual .mmi files for QCJob.
For NuQCJob, minimap2_databases is expected to be the path to a directory containing two subdirectories: 'metagenomic'
and 'metatranscriptomic'. Each directory should contain or symlink to the appropriate .mmi files needed for that Assay
type.
