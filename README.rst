# Sequence Processing Pipeline

A Jupyter notebook to assist wet lab shotgun pipeline.
A packaged Python-based implementation of Knight Lab's sequencing process pipeline.

## Installation

To install this package, first clone these repositories from GitHub:

```bash
git clone https://github.com/biocore/mg-scripts.git
git clone https://github.com/biocore/metagenomics_pooling_notebook.git
```

Create a Python3 Conda environment in which to run the notebook:

```bash
conda create -n sp_pipeline 'python==3.9' numpy pandas click scipy matplotlib 
```

Activate the Conda environment:

```bash
source activate sp_pipeline
```

Change directory to the downloaded repository folder and install:

```bash
cd metagenomics_pooling_notebook
pip install -e '.[all]'
```

Change directory to the parent folder:

```bash
cd ..
```

Change directory to the downloaded repository folder and install:

```bash
cd mg-scripts
pip install -e .
```

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
