from os.path import join, realpath, dirname
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.PipelineError import PipelineError
from functools import partial
import unittest


class TestConvertJob(unittest.TestCase):
    def test_creation(self):
        # currently, the root directory of a job is going to be the
        # root directory of the single sequence directory given as
        # input by the user. Later we can re-introduce directories
        # that contain multiple BCL root directories.
        path = partial(join, 'sequence_processing_pipeline', 'tests', 'data')
        run_dir = path('sample-sequence-directory')
        sample_sheet_path = path('good-sample-sheet.csv')
        output_directory = path('output_directory')
        inv_input_directory = path('inv_input_directory')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        # ConvertJob should assert due to invalid_input_directory.
        with self.assertRaises(PipelineError):
            ConvertJob(inv_input_directory, sample_sheet_path,
                       output_directory, 'qiita', 1, 16, 24, '10gb',
                       'tests/bin/bcl-convert', [], qiita_id)

        job = ConvertJob(run_dir, sample_sheet_path, output_directory, 'qiita',
                         1, 16, 24, '10gb', 'tests/bin/bcl-convert', [],
                         qiita_id)
        with open(job._generate_job_script()) as f:
            obs = ''.join(f.readlines())

        full_run_dir = join(
            realpath(dirname(__file__)), 'data', 'sample-sequence-directory')

        self.assertEqual(obs, SCRIPT_EXP.format(full_run_dir=full_run_dir))


SCRIPT_EXP = ''.join([
    '#!/bin/bash\n',
    '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_ConvertJob\n',
    '#PBS -q qiita\n',
    '#PBS -l nodes=1:ppn=16\n',
    '#PBS -V\n',
    '#PBS -l walltime=24:00:00\n',
    '#PBS -m bea\n',
    '#PBS -M ccowart@ucsd.edu,qiita.help@gmail.com\n',
    '#PBS -l pmem=10gb\n',
    '#PBS -o localhost:{full_run_dir}/ConvertJob.out.log\n',
    '#PBS -e localhost:{full_run_dir}/ConvertJob.err.log\n',
    'set -x\n',
    'cd {full_run_dir}\n',
    'module load \n',
    'tests/bin/bcl-convert --sample-sheet sequence_processing_pipeline/tests/'
    'data/good-sample-sheet.csv --output-directory '
    'sequence_processing_pipeline/tests/data/output_directory '
    '--bcl-input-directory . --bcl-num-decompression-threads 8 '
    '--bcl-num-conversion-threads 8 --bcl-num-compression-threads 16 '
    '--bcl-num-parallel-tiles 8 --bcl-sampleproject-subdirectories '
    'true --force\n'])

if __name__ == '__main__':
    unittest.main()
