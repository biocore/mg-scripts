import json
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import unittest
from os import utime, makedirs
from os.path import join
from copy import deepcopy
from time import time


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
        self.good_output_file_path = ('sequence_processing_pipeline/'
                                     'my_output_file_path')

    def test_creation(self):
        # Pipeline should assert due to config_file
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.bad_config_file,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        self.assertEqual(str(e.exception), "'search_paths' is not a key in "
                                           "sequence_processing_pipeline/tests"
                                           "/data/bad_configuration.json")

        # Pipeline should assert due to an invalid config file path.
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.invalid_config_file,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        self.assertEqual(str(e.exception), 'does/not/exist/configuration.json '
                                           'does not exist.')

        # Pipeline should assert on config_file = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(None,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        self.assertEqual(str(e.exception), 'configuration_file_path cannot be '
                                           'None')

        # Pipeline should assert due to invalid_run_id
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     self.invalid_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        self.assertEqual(str(e.exception), "A run-dir for 'not-sample-sequence"
                                           "-directory' could not be found")

        # Pipeline should assert on run_id = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     None,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        with open(join('sequence_processing_pipeline', 'configuration.json'),
                  'r') as f:
            cfg = json.load(f)
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = -1
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_output_file_path,
                         'my_qiita_id',
                         cpy_cfg)

            self.assertEqual(str(e.exception), 'older_than and younger_than '
                                               'cannot be less than zero.')
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['older_than'] = -1
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_output_file_path,
                         'my_qiita_id',
                         cpy_cfg)
            self.assertEqual(str(e.exception), 'older_than and younger_than '
                                               'cannot be less than zero.')
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = 20
                cpy_cfg['configuration']['pipeline']['older_than'] = 30
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_output_file_path,
                         'my_qiita_id',
                         cpy_cfg)
            self.assertEqual(str(e.exception), 'older_than cannot be equal to '
                                               'or less than younger_than.')

    def test_filter_directories_for_time(self):
        # get a good configuration dict from file.
        with open(join('sequence_processing_pipeline', 'configuration.json'),
                  'r') as f:
            cfg = json.load(f)
            # set a range that sample-sequence-directory's timestamps must
            # fall within: between 1 and 2 hours.
            cfg['configuration']['pipeline']['younger_than'] = 2
            cfg['configuration']['pipeline']['older_than'] = 1

            # Set directory's timestamp to between one and two hours and
            # verify it is returned by filter_directories_for_time().

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value older than 1 hour ago + 5 min.
            older_than = current_time - (3600 + (5 * 60))
            tp = join(self.data_directory, 'sample-sequence-directory')
            utime(tp, (older_than, older_than))

            pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                self.good_output_file_path, 'my_qiita_id',
                                cfg)

            obs = pipeline.is_within_time_range(tp)
            self.assertEqual(obs, True)

            # Set directory's timestamp to just older than two hours and
            # verify it is not returned by filter_directories_for_time().

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value older than 2 hour ago + 5 min.
            older_than = current_time - (7200 + (5 * 60))
            utime(tp, (older_than, older_than))

            obs = pipeline.is_within_time_range(tp)
            self.assertEqual(obs, False)

            # Set directory's timestamp to just under one hour and
            # verify it is not returned by filter_directories_for_time().

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value younger than 1 hour ago.
            older_than = current_time - 3300
            utime(join(self.data_directory, 'sample-sequence-directory'),
                  (older_than, older_than))

            obs = pipeline.is_within_time_range(tp)
            self.assertEqual(obs, False)


if __name__ == '__main__':
    unittest.main()
