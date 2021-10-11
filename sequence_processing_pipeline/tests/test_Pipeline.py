from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import logging
import unittest
from os import makedirs


logging.basicConfig(level=logging.DEBUG)


class TestPipeline(unittest.TestCase):
    def test_creation(self):
        makedirs('sequence_processing_pipeline/tests/data/sequencing'
                 '/knight_lab_completed_runs', exist_ok=True)

        configuration_file = 'sequence_processing_pipeline/configuration.json'

        # shutil.rmtree(
        # 'sequence_processing_pipeline/tests/data/invalid_output_directory')
        # shutil.rmtree(('sequence_processing_pipeline/tests'
        #               '/data/invalid_output_directory_two'),
        #              ignore_errors=False)

        # test invalid input_directory
        self.assertRaises(PipelineError,
                          Pipeline,
                          configuration_file,
                          'not-valid-sample-sequence-directory')

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
        # self.assertRaises(PipelineError, Pipeline, input_directory,
        #                  output_directory, final_output_directory, nprocs=18)

        # nprocs must be greater than 0.
        # self.assertRaises(PipelineError, Pipeline, input_directory,
        #                  output_directory, final_output_directory, nprocs=0)

        # window of time must be legitimate.
        # (older_than can't be less than or equal to younger_than.)
        # self.assertRaises(PipelineError, Pipeline, input_directory,
        #                   output_directory, final_output_directory,
        #                   younger_than=24, older_than=24)

        # using all default parameters should not raise an Error.
        # test_path = partial(join, dirname(abspath(__file__)))
        # input_directory = test_path('data/sample-sequence-directory')

        # msg = None
        # try:
        #     Pipeline(configuration_file, input_directory, run_id='run-id')
        # except PipelineError as e:
        #     msg = str(e)

        # self.assertEqual(msg, None)
