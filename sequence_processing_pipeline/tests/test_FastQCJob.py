import unittest
from os.path import join, exists
from functools import partial
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.PipelineError import PipelineError
from os import makedirs
import shutil


class TestFastQCJob(unittest.TestCase):
    def setUp(self):
        package_root = 'sequence_processing_pipeline'
        self.path = partial(join, package_root, 'tests', 'data')
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.maxDiff = None
        self.output_path = self.path('output_dir')
        self.raw_fastq_files_path = ('sequence_processing_pipeline/tests/data'
                                     '/211021_A00000_0000_SAMPLE/Data/Fastq/p'
                                     'roject1')
        self.processed_fastq_files_path = ('sequence_processing_pipeline/tests'
                                           '/data/211021_A00000_0000_SAMPLE/sa'
                                           'mple-sequence-directory')
        self.config_yml = join(package_root, 'multiqc-bclconvert-config.yaml')
        self.qc_root_path = join(self.output_path, 'QCJob')
        makedirs(self.qc_root_path, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self.qc_root_path)

    def test_config_file_not_found(self):
        with self.assertRaises(PipelineError) as e:
            FastQCJob(self.qc_root_path, self.output_path,
                      self.raw_fastq_files_path,
                      self.processed_fastq_files_path, 16, 16,
                      'sequence_processing_pipeline/tests/bin/fastqc', [],
                      self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                      ('sequence_processing_pipeline/'
                       'not-multiqc-bclconvert-config.yaml'), 1000)

        self.assertEqual(str(e.exception), "file 'sequence_processing_pipeline"
                                           "/not-multiqc-bclconvert-config."
                                           "yaml' does not exist.")

    def test_generate_job_scripts(self):
        job = FastQCJob(self.qc_root_path, self.output_path,
                        self.raw_fastq_files_path,
                        self.processed_fastq_files_path,
                        16, 16,
                        'sequence_processing_pipeline/tests/bin/fastqc', [],
                        self.qiita_job_id, 'queue_name', 4, 23, '8g', 30,
                        self.config_yml, 1000)

        self.assertEqual(exists(join(job.output_path, 'FastQCJob.sh')), True)
        self.assertEqual(exists(join(job.output_path,
                                     'FastQCJob.array-details')), True)


if __name__ == '__main__':
    unittest.main()
