import logging
from sequence_processing_pipeline.BCL2FASTQJob import BCL2FASTQJob
from os.path import join, abspath
import re
from os import makedirs
import unittest
from sequence_processing_pipeline.Pipeline import Pipeline
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import exists


logging.basicConfig(level=logging.DEBUG)



class TestBCL2FASTQJob(unittest.TestCase):
    def test_creation(self):
        # currently, the root directory of a job is going to be the
        # root directory of the single sequence directory given as
        # input by the user. Later we can re-introduce directories
        # that contain multiple BCL root directories.
        root_dir = 'tests/data/sample-sequence-directory'
        sample_sheet_path = 'tests/data/good-sample-sheet.csv'
        output_directory = 'tests/data/output_directory'
        bcl2fastq_path = 'tests/data/fake_bcl2fast2.sh'
        invalid_input_directory = 'tests/data/invalid_input_directory'
        invalid_output_directory = 'tests/data/invalid_output_directory'
        invalid_final_output_directory = 'tests/data/invalid_output_directory_two'

        # BCL2FASTQJob should assert due to invalid_input_directory.
        self.assertRaises(PipelineError, BCL2FASTQJob, invalid_input_directory, sample_sheet_path, output_directory, bcl2fastq_path)
        # BCL2FASTQJob should assert due root_dir/Data directory being devoid of BCL files.
        self.assertRaises(PipelineError, BCL2FASTQJob, root_dir, sample_sheet_path, output_directory, bcl2fastq_path)
        # Create fake BCL files in the directory structure seen in real-world examples.
        makedirs('tests/data/sample-sequence-directory/Data/Intensities/BaseCalls/L003', exist_ok=True)
        with open('tests/data/sample-sequence-directory/Data/Intensities/BaseCalls/L003/fake.bcl', 'w')  as f:
            f.write('this is a text file.')

        # BCL2FASTQJob should not assert an error now.
        msg = None
        try:
            BCL2FASTQJob(root_dir, sample_sheet_path, output_directory, bcl2fastq_path)
        except PipelineError as e:
            msg = str(e)

        self.assertEqual(msg, None)

        '''
        sample_sheet_path is not a real sheet path
        sample_sheet path is not a valid sheet
        (we can assume has the data points we need if it passes)
        output directory must be creatable if not present.
        bcl2fastq_path needs to point to a real executable binary.
        '''


