#!/usr/bin/env python
from setuptools import setup
from glob import glob


__version__ = "2021.08.24"


classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""

with open('README.md') as f:
    long_description = f.read()

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='sequence-processing-pipeline',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Sequence processing pipeline',
      author="Qiita development team",
      author_email="qiita.help@gmail.com",
      url='https://github.com/biocore/mg-scripts',
      packages=['sequence_processing_pipeline'],
      scripts=glob('scripts/*'),
      setup_requires=['nose3', 'numpy', 'cython'],
      install_requires=[
        'click', 'requests', 'pandas', 'flake8',
        'metapool @ https://github.com/charles-cowart/metagenomics_pooling_notebook-1/'
        'archive/master.zip'],
      classifiers=classifiers
      )
