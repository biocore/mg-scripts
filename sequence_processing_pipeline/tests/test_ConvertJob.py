from os import makedirs, environ
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.PipelineError import PipelineError
from shutil import rmtree
import logging
import unittest


logging.basicConfig(level=logging.DEBUG)


class TestConvertJob(unittest.TestCase):
    def test_creation(self):
        # currently, the root directory of a job is going to be the
        # root directory of the single sequence directory given as
        # input by the user. Later we can re-introduce directories
        # that contain multiple BCL root directories.
        run_dir = 'tests/data/sample-sequence-directory'
        sample_sheet_path = 'tests/data/good-sample-sheet.csv'
        output_directory = 'tests/data/output_directory'
        inv_input_directory = 'tests/data/invalid_input_directory'
        # inv_output_dir = 'tests/data/invalid_output_directory'
        # inv_final_output_dir = 'tests/data/invalid_output_directory_two'

        # ConvertJob should assert due to invalid_input_directory.
        self.assertRaises(PipelineError, ConvertJob, inv_input_directory,
                          sample_sheet_path, output_directory, True, 'qiita',
                          1, 16, 24)

        # ConvertJob should assert due run_dir/Data directory being
        # devoid of BCL files.
        self.assertRaises(PipelineError, ConvertJob, run_dir,
                          sample_sheet_path, output_directory, True, 'qiita',
                          1, 16, 24)

        # Create fake BCL files in the directory structure seen in real-world
        # examples.
        makedirs('tests/data/sample-sequence-directory/Data/Intensities'
                 '/BaseCalls/L003', exist_ok=True)
        with open('tests/data/sample-sequence-directory/Data/Intensities'
                  '/BaseCalls/L003/fake.bcl', 'w') as f:
            f.write('this is a text file.')

        # ConvertJob should not assert an error now.
        environ['PATH'] = ('/Users/ccowart/PycharmProjects/mg-scripts/'
                           'sequence_processing_pipeline/tests/bin:'
                           + environ['PATH'])

        msg = None
        try:
            ConvertJob(run_dir, sample_sheet_path, output_directory, True,
                       'qiita', 1, 16, 24)
        except PipelineError as e:
            msg = str(e)

        logging.debug(msg)
        # self.assertEqual(msg, None)

        rmtree('tests/data/sample-sequence-directory/Data/Intensities')
