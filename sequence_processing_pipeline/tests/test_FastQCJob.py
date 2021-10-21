import unittest
from os.path import join
from functools import partial
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.PipelineError import PipelineError


class TestFastQCJob(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.path = partial(join, 'sequence_processing_pipeline', 'tests',
                            'data')
        self.run_dir = self.path('sample-sequence-directory')
        self.output_directory = self.path('output_directory')
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

    def test_fastqc_not_found(self):
        with self.assertRaises(PipelineError) as e:
            FastQCJob(self.run_dir, self.output_directory, 16, 16,
                      'sequence_processing_pipeline/tests/bin/not-fastqc',
                      [], self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                      ('sequence_processing_pipeline/'
                       'multiqc-bclconvert-config.yaml'))

        self.assertEqual(str(e.exception), "file 'sequence_processing_pipeline"
                                           "/tests/bin/not-fastqc' does not ex"
                                           "ist.")

    def test_config_file_not_found(self):
        with self.assertRaises(PipelineError) as e:
            FastQCJob(self.run_dir, self.output_directory, 16, 16,
                      'sequence_processing_pipeline/tests/bin/fastqc',
                      [], self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                      ('sequence_processing_pipeline/'
                       'not-multiqc-bclconvert-config.yaml'))

        self.assertEqual(str(e.exception), "file 'sequence_processing_pipeline"
                                           "/not-multiqc-bclconvert-config."
                                           "yaml' does not exist.")

    def test_generate_job_script_path(self):
        job = FastQCJob(self.run_dir, self.output_directory, 16, 16,
                        'sequence_processing_pipeline/tests/bin/fastqc',
                        [], self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                        'sequence_processing_pipeline/'
                        'multiqc-bclconvert-config.yaml')

        exp = (f'{self.run_dir}/FastQCJob_2.sh',
               f'localhost:/{self.run_dir}/FastQCJob_2.out.log',
               f'localhost:/{self.run_dir}/FastQCJob_2.err.log')

        self.assertEqual(job.generate_job_script_path(), exp)


if __name__ == '__main__':
    unittest.main()
