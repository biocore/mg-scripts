import unittest
from sequence_processing_pipeline.FastQCJob import FastQCJob
from os.path import join
from functools import partial


class TestFastQCJob(unittest.TestCase):
    def test_creation(self):
        path = partial(join, 'sequence_processing_pipeline', 'tests', 'data')
        run_dir = path('sample-sequence-directory')
        output_directory = path('output_directory')
        qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        with self.assertRaises(FileNotFoundError):
            FastQCJob(run_dir, output_directory, '8', '16',
                      '/path/fastqc_path', [], qiita_job_id, 'not_run_id',
                      'queue_name', '4', '23:00:00', '8g')

        job = FastQCJob(run_dir, output_directory, '8', '16',
                        '/path/fastqc_path', [], qiita_job_id, 'run-id',
                        'queue_name', '4', '23:00:00', '8g')
        exp = (
            f'{run_dir}/FastQCJob_1.sh',
            f'localhost:/{run_dir}/FastQCJob_1.out.log',
            f'localhost:/{run_dir}/FastQCJob_1.err.log')

        self.assertEqual(job.generate_job_script_path(), exp)


if __name__ == '__main__':
    unittest.main()
