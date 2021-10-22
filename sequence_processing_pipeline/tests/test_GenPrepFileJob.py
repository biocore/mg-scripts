from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from os.path import join
from functools import partial
import unittest


class TestGenPrepFileJob(unittest.TestCase):
    def test_creation(self):
        package_root = 'sequence_processing_pipeline'
        self.path = partial(join, package_root, 'tests', 'data')
        run_dir = self.path('sample-sequence-directory')
        sample_sheet_path = self.path('good-sample-sheet.csv')
        self.output_path = self.path('output_dir')
        self.maxDiff = None

        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        job = GenPrepFileJob(run_dir, self.output_path, sample_sheet_path,
                             'seqpro', [], qiita_id)
        exp = ['seqpro', run_dir, sample_sheet_path,
               join(self.output_path, 'GenPrepFileJob')]

        self.assertEqual(job.command, exp)


if __name__ == '__main__':
    unittest.main()
