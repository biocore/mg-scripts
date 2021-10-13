import unittest
from os.path import join
from functools import partial
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.PipelineError import PipelineError

import logging
logging.basicConfig(level=logging.DEBUG)


class TestFastQCJob(unittest.TestCase):
    def test_creation(self):
        path = partial(join, 'sequence_processing_pipeline', 'tests', 'data')
        run_dir = path('sample-sequence-directory')
        output_directory = path('output_directory')
        qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        with self.assertRaises(PipelineError) as e:
            FastQCJob(run_dir, output_directory, 16, 16,
                      'sequence_processing_pipeline/tests/bin/not-fastqc',
                      [], qiita_job_id, 'queue_name', 4, 23, '8g', 30)

        self.assertEqual(str(e.exception), "file 'sequence_processing_pipeline"
                                           "/tests/bin/not-fastqc' does not ex"
                                           "ist.")

        job = FastQCJob(run_dir, output_directory, 16, 16,
                        'sequence_processing_pipeline/tests/bin/fastqc',
                        [], qiita_job_id, 'queue_name', 4, 23, '8g', 30)

        exp = (f'{run_dir}/FastQCJob_2.sh',
               f'localhost:/{run_dir}/FastQCJob_2.out.log',
               f'localhost:/{run_dir}/FastQCJob_2.err.log')

        self.assertEqual(job.generate_job_script_path(), exp)


if __name__ == '__main__':
    unittest.main()
