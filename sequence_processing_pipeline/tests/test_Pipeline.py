import shutil

from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import logging
import unittest
from os.path import exists


logging.basicConfig(level=logging.DEBUG)


class TestPipeline(unittest.TestCase):
    def test_creation(self):
        input_directory = 'tests/data/sample-sequence-directory'
        output_directory = 'tests/data/output_directory'
        final_output_directory = 'tests/data/final_output_directory'
        inv_input_directory = 'tests/data/invalid_input_directory'
        inv_output_directory = 'tests/data/invalid_output_directory'
        inv_final_output_directory = 'tests/data/invalid_output_directory_two'

        # all parameters defined
        # pipeline = Pipeline(input_directory,
        #                    output_directory,
        #                    final_output_directory,
        #                    younger_than=48,
        #                    older_than=24,
        #                    nprocs=16)

        # test output_directory == final_output_directory
        self.assertRaises(PipelineError,
                          Pipeline,
                          input_directory,
                          output_directory,
                          output_directory)

        # shutil.rmtree('sequence_processing_pipeline/tests/data/invalid_output_directory')
        shutil.rmtree(
            'sequence_processing_pipeline/tests/data/invalid_output_directory_two')
        exit(0)
        # test invalid input_directory
        self.assertRaises(PipelineError,
                          Pipeline,
                          inv_input_directory,
                          output_directory, final_output_directory)

        # test non-existant output_directory
        # pipeline = Pipeline(input_directory,
        #                    inv_output_directory,
        #                    final_output_directory)
        # self.assertTrue(exists(inv_output_directory))

        # test non-existant final_output_directory
        # pipeline = Pipeline(input_directory,
        #                    output_directory,
        #                    inv_final_output_directory)
        # self.assertTrue(exists(inv_final_output_directory))

        # nprocs exceeds 16 nproc limit
        self.assertRaises(PipelineError, Pipeline, input_directory,
                          output_directory, final_output_directory, nprocs=18)

        # nprocs must be greater than 0.
        self.assertRaises(PipelineError, Pipeline, input_directory,
                          output_directory, final_output_directory, nprocs=0)

        # window of time must be legitimate.
        # (older_than can't be less than or equal to younger_than.)
        self.assertRaises(PipelineError, Pipeline, input_directory,
                          output_directory, final_output_directory,
                          younger_than=24, older_than=24)

        # using all default parameters should not raise an Error.
        msg = None
        try:
            Pipeline(input_directory, output_directory, final_output_directory)
        except PipelineError as e:
            msg = str(e)

        self.assertEqual(msg, None)
