from os import makedirs
from sequence_processing_pipeline.BCLConvertJob import BCLConvertJob
from sequence_processing_pipeline.PipelineError import PipelineError
from shutil import rmtree
import logging
import unittest


logging.basicConfig(level=logging.DEBUG)


class TestBCLConvertJob(unittest.TestCase):
    def test_creation(self):
        # currently, the root directory of a job is going to be the
        # root directory of the single sequence directory given as
        # input by the user. Later we can re-introduce directories
        # that contain multiple BCL root directories.
        root_dir = 'tests/data/sample-sequence-directory'
        sample_sheet_path = 'tests/data/good-sample-sheet.csv'
        output_directory = 'tests/data/output_directory'
        bcl2fastq_path = 'tests/data/fake_bcl2fastq'
        invalid_input_directory = 'tests/data/invalid_input_directory'
        invalid_output_directory = 'tests/data/invalid_output_directory'
        invalid_final_output_directory = 'tests/data/invalid_output_directory_two'

        # ConvertConvertBCL2FastqJob should assert due to invalid_input_directory.
        self.assertRaises(PipelineError, BCLConvertJob, invalid_input_directory, sample_sheet_path, output_directory, bcl2fastq_path)

        # ConvertConvertBCL2FastqJob should assert due root_dir/Data directory being devoid of BCL files.
        self.assertRaises(PipelineError, BCLConvertJob, root_dir, sample_sheet_path, output_directory, bcl2fastq_path)

        # Create fake BCL files in the directory structure seen in real-world examples.
        makedirs('tests/data/sample-sequence-directory/Data/Intensities/BaseCalls/L003', exist_ok=True)
        with open('tests/data/sample-sequence-directory/Data/Intensities/BaseCalls/L003/fake.bcl', 'w')  as f:
            f.write('this is a text file.')

        # ConvertConvertBCL2FastqJob should not assert an error now.
        msg = None
        try:
            BCLConvertJob(root_dir, sample_sheet_path, output_directory, bcl2fastq_path)
        except PipelineError as e:
            msg = str(e)

        logging.debug(msg)
        self.assertEqual(msg, None)

        rmtree('tests/data/sample-sequence-directory/Data/Intensities')

