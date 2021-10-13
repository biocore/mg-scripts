import json
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import logging
import unittest
from os import utime, makedirs
from os.path import join
from copy import deepcopy
from time import time


logging.basicConfig(level=logging.DEBUG)


class TestPipeline(unittest.TestCase):
    def setUp(self):
        makedirs('sequence_processing_pipeline/tests/data/sequencing'
                 '/knight_lab_completed_runs', exist_ok=True)

        self.good_config_file = ('sequence_processing_pipeline/'
                                 'configuration.json')
        self.bad_config_file = ('sequence_processing_pipeline/tests/data/'
                                'bad_configuration.json')
        self.invalid_config_file = 'does/not/exist/configuration.json'
        self.good_run_id = 'sample-sequence-directory'
        self.invalid_run_id = 'not-sample-sequence-directory'
        self.data_directory = 'sequence_processing_pipeline/tests/data'
        self.maxDiff = None

    def test_creation(self):
        # Pipeline should assert due to config_file
        with self.assertRaises(PipelineError):
            Pipeline(self.bad_config_file, self.good_run_id)

        # Pipeline should assert due to an invalid config file path.
        with self.assertRaises(PipelineError):
            Pipeline(self.invalid_config_file, self.good_run_id)

        # Pipeline should assert on config_file = None
        with self.assertRaises(PipelineError):
            Pipeline(None, self.good_run_id)

        # Pipeline should assert due to invalid_run_id
        with self.assertRaises(PipelineError):
            Pipeline(self.good_config_file, self.invalid_run_id)

        # Pipeline should assert on run_id = None
        with self.assertRaises(PipelineError):
            Pipeline(self.good_config_file, None)

        with open(join('sequence_processing_pipeline', 'configuration.json'),
                  'r') as f:
            cfg = json.load(f)
            with self.assertRaises(PipelineError):
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = -1
                logging.debug(json.dumps(cpy_cfg, indent=2))
                Pipeline(self.good_config_file, self.good_run_id,
                         config_dict=cpy_cfg)
            with self.assertRaises(PipelineError):
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['older_than'] = -1
                Pipeline(self.good_config_file, self.good_run_id,
                         config_dict=cpy_cfg)
            with self.assertRaises(PipelineError):
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = 20
                cpy_cfg['configuration']['pipeline']['older_than'] = 30
                Pipeline(self.good_config_file, self.good_run_id,
                         config_dict=cpy_cfg)

    def test_find_bcl_directories(self):
        pipeline = Pipeline(self.good_config_file, self.good_run_id)
        obs = pipeline.find_bcl_directories(self.data_directory)
        # there are a number of directories under 'data'. Only
        # 'sample-sequence-directory' has the layout expected of a
        # run-directory, however.
        exp = ["sequence_processing_pipeline/tests/data/"
               "sample-sequence-directory"]
        self.assertEqual(obs, exp)

    def test_filter_directories_for_time(self):
        # get a good configuration dict from file.
        with open(join('sequence_processing_pipeline', 'configuration.json'),
                  'r') as f:
            cfg = json.load(f)
            # set a range that sample-sequence-directory's timestamps must
            # fall within: between 1 and 2 hours.
            cfg['configuration']['pipeline']['younger_than'] = 2
            cfg['configuration']['pipeline']['older_than'] = 1

            # Test 1

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value older than 1 hour ago + 5 min.
            older_than = current_time - (3600 + (5 * 60))
            utime(join(self.data_directory, 'sample-sequence-directory'),
                  (older_than, older_than))

            pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                config_dict=cfg)
            # there should be only sample-sequence-directory in the results.
            results = pipeline.find_bcl_directories(self.data_directory)

            obs = pipeline.filter_directories_for_time(results)
            exp = ["sequence_processing_pipeline/tests/data/"
                   "sample-sequence-directory"]
            self.assertEqual(obs, exp)

            # Test 2

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value older than 2 hour ago + 5 min.
            older_than = current_time - (7200 + (5 * 60))
            utime(join(self.data_directory, 'sample-sequence-directory'),
                  (older_than, older_than))

            pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                config_dict=cfg)
            # there should be only sample-sequence-directory in the results.
            results = pipeline.find_bcl_directories(self.data_directory)

            obs = pipeline.filter_directories_for_time(results)
            exp = []
            self.assertEqual(obs, exp)

            # Test 3

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value younger than 1 hour ago.
            older_than = current_time - 3300
            utime(join(self.data_directory, 'sample-sequence-directory'),
                  (older_than, older_than))

            pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                config_dict=cfg)
            # there should be only sample-sequence-directory in the results.
            results = pipeline.find_bcl_directories(self.data_directory)

            obs = pipeline.filter_directories_for_time(results)
            exp = []
            self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
