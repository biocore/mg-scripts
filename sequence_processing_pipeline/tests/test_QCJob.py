import unittest
from os.path import join
from functools import partial
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.PipelineError import PipelineError


import logging
logging.basicConfig(level=logging.DEBUG)


class TestQCJob(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.path = partial(join, 'sequence_processing_pipeline', 'tests',
                            'data')
        self.run_dir = self.path('good-sample-sheet.csv')
        self.sample_sheet_path = ''
        self.mmi_db_path = self.path('mmi.db')
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.products_dir = self.path('sample-sequence-directory',
                                      'sample-sequence-directory')

    def test_qcjob_creation(self):
        with self.assertRaises(PipelineError) as e:
            QCJob(self.run_dir, 'not/path/to/sample/sheet', self.mmi_db_path,
                  'queue_name', 1, 16, 24, '8gb', 'fastp', 'minimap2',
                  'samtools', [], self.qiita_job_id, 30, self.products_dir)

        self.assertEqual(str(e.exception), "file 'not/path/to/sample/sheet' "
                                           "does not exist.")


if __name__ == '__main__':
    unittest.main()
