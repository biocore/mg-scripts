from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.PipelineError import PipelineError
import logging
import unittest


logging.basicConfig(level=logging.DEBUG)


class TestGenPrepFileJob(unittest.TestCase):
    def test_creation(self):
        run_dir = 'tests/data/sample-sequence-directory'
        sample_sheet_path = 'tests/data/good-sample-sheet.csv'
        output_directory = 'tests/data/output_directory'
        qiita_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        msg = None
        try:
            GenPrepFileJob(run_dir, sample_sheet_path, output_directory,
                           'seqpro', [], qiita_id)
        except PipelineError as e:
            msg = str(e)

        logging.debug(msg)
        # self.assertEqual(msg, None)
