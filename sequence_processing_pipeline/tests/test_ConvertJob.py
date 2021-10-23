from os.path import join
from os import makedirs
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
        inv_input_directory = path('inv_input_directory')
        good_output_path = path('output_dir')
        makedirs(good_output_path, exist_ok=True)
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.maxDiff = None

        # ConvertJob should assert due to invalid_input_directory.
        with self.assertRaises(PipelineError) as e:
            ConvertJob(inv_input_directory, good_output_path,
                       sample_sheet_path, 'qiita', 1, 16, 24, '10gb',
                       'tests/bin/bcl-convert', [], qiita_id)

        self.assertEqual(str(e.exception), ("directory_path '"
                                            f"{path('inv_input_directory')}' "
                                            "does not exist."))

        job = ConvertJob(run_dir, good_output_path, sample_sheet_path,
                         'qiita', 1, 16, 24, '10gb', 'tests/bin/bcl-convert',
                         [], qiita_id)

        job._generate_job_script()

        with open(join(good_output_path, 'ConvertJob', 'ConvertJob.sh')) as f:
            obs = ''.join(f.readlines())

        self.assertEqual(obs, SCRIPT_EXP.format(gop=good_output_path,
                                                run_dir=run_dir))


SCRIPT_EXP = ''.join([
    '#!/bin/bash\n',
    '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_ConvertJob\n',
    '#PBS -q qiita\n',
    '#PBS -l nodes=1:ppn=16\n',
    '#PBS -V\n',
    '#PBS -l walltime=24:00:00\n',
    '#PBS -m bea\n',
    '#PBS -M qiita.help@gmail.com\n',
    '#PBS -l pmem=10gb\n',
    '#PBS -o localhost:{gop}/ConvertJob/logs/qsub_stdout.log\n',
    '#PBS -e localhost:{gop}/ConvertJob/logs/qsub_stderr.log\n',
    'set -x\n',
    'date\n',
    'hostname\n',
    'cd {run_dir}\n',
    'tests/bin/bcl-convert --sample-sheet sequence_processing_pipeline/tests/'
    'data/good-sample-sheet.csv --output-directory {gop}/ConvertJob '
    '--bcl-input-directory . '
    '--bcl-num-decompression-threads 16 --bcl-num-conversion-threads 16 '
    '--bcl-num-compression-threads 16 --bcl-num-parallel-tiles 16 '
    '--bcl-sampleproject-subdirectories true --force\n'])


if __name__ == '__main__':
    unittest.main()
