from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from os.path import join
from functools import partial
import unittest


class TestGenPrepFileJob(unittest.TestCase):
    def test_creation(self):
        path = partial(join, 'sequence_processing_pipeline', 'tests', 'data')
        run_dir = path('sample-sequence-directory')
        sample_sheet_path = path('good-sample-sheet.csv')
        output_directory = path('output_directory')
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        job = GenPrepFileJob(run_dir, sample_sheet_path, output_directory,
                             'seqpro', [], qiita_id)
        exp = ['seqpro', run_dir, sample_sheet_path,
               join(output_directory, 'prep_files')]

        self.assertEqual(job.command, exp)


if __name__ == '__main__':
    unittest.main()
