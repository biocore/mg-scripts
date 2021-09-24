from os import makedirs
from sequence_processing_pipeline.ConvertBCL2FastqJob import ConvBCL2FastqJob
from sequence_processing_pipeline.PipelineError import PipelineError
from shutil import rmtree
import logging
import unittest


logging.basicConfig(level=logging.DEBUG)


class TestConvertBCL2FastqJob(unittest.TestCase):
    def test_creation(self):
        # currently, the root directory of a job is going to be the
        # root directory of the single sequence directory given as
        # input by the user. Later we can re-introduce directories
        # that contain multiple BCL root directories.
        run_dir = 'tests/data/sample-sequence-directory'
        sample_sheet_path = 'tests/data/good-sample-sheet.csv'
        output_directory = 'tests/data/output_directory'
        bcl2fastq_path = 'tests/data/fake_bcl2fastq'
        inv_input_directory = 'tests/data/invalid_input_directory'
        # inv_output_dir = 'tests/data/invalid_output_directory'
        # inv_final_output_dir = 'tests/data/invalid_output_directory_two'

        # ConvertBCL2FastqJob should assert due to invalid_input_directory.
        self.assertRaises(PipelineError, ConvBCL2FastqJob, inv_input_directory,
                          sample_sheet_path, output_directory, bcl2fastq_path)

        # ConvertBCL2FastqJob should assert due run_dir/Data directory being
        # devoid of BCL files.
        self.assertRaises(PipelineError, ConvBCL2FastqJob, run_dir,
                          sample_sheet_path, output_directory, bcl2fastq_path)

        # Create fake BCL files in the directory structure seen in real-world
        # examples.
        makedirs('tests/data/sample-sequence-directory/Data/Intensities'
                 '/BaseCalls/L003', exist_ok=True)
        with open('tests/data/sample-sequence-directory/Data/Intensities'
                  '/BaseCalls/L003/fake.bcl', 'w') as f:
            f.write('this is a text file.')

        # ConvertBCL2FastqJob should not assert an error now.
        msg = None
        try:
            ConvBCL2FastqJob(run_dir, sample_sheet_path, output_directory,
                             bcl2fastq_path)
        except PipelineError as e:
            msg = str(e)

        logging.debug(msg)
        self.assertEqual(msg, None)

        rmtree('tests/data/sample-sequence-directory/Data/Intensities')
