import json
import os

from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import unittest
from os import utime, makedirs
from os.path import basename, join, abspath
from copy import deepcopy
from time import time
from functools import partial


class TestPipeline(unittest.TestCase):
    def setUp(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.bad_config_file = self.path('bad_configuration.json')
        self.invalid_config_file = 'does/not/exist/configuration.json'
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.invalid_run_id = 'not-sample-sequence-directory'
        self.good_output_file_path = self.path('output_dir')
        makedirs(self.good_output_file_path, exist_ok=True)
        self.maxDiff = None
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.bad_sample_sheet_path = self.path('duplicate_sample-sample-sheet'
                                               '.csv')
        self.runinfo_file = self.path(self.good_run_id, 'RunInfo.xml')
        self.rtacomplete_file = self.path(self.good_run_id, 'RTAComplete.txt')

        # most of the tests here were written with the assumption that these
        # files already exist.
        self.create_runinfo_file()
        self.create_rtacomplete_file()

    def tearDown(self):
        # Pipeline is now the only class aware of these files, hence they
        # can be deleted at the end of testing.
        self.delete_runinfo_file()
        self.delete_rtacomplete_file()

    def create_runinfo_file(self):
        with open(self.runinfo_file, 'w') as f:
            f.write("")

    def delete_runinfo_file(self):
        try:
            os.remove(self.runinfo_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def create_rtacomplete_file(self):
        with open(self.rtacomplete_file, 'w') as f:
            f.write("")

    def delete_rtacomplete_file(self):
        try:
            os.remove(self.rtacomplete_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def test_required_file_checks(self):
        # begin this test by deleting the RunInfo.txt file and verifying that
        # Pipeline object will raise an Error.
        self.delete_runinfo_file()
        with self.assertRaisesRegex(PipelineError, "required file 'RunInfo.xml"
                                                   "' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_output_file_path, 'my_qiita_id', None)

        # delete RTAComplete.txt and recreate RunInfo.txt file to verify that
        # an Error is raised when only RTAComplete.txt is missing.
        self.delete_rtacomplete_file()
        self.create_runinfo_file()
        with self.assertRaisesRegex(PipelineError, "required file 'RTAComplete"
                                                   ".txt' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_output_file_path, 'my_qiita_id', None)

    def test_creation(self):
        # Pipeline should assert due to config_file
        with self.assertRaisesRegex(PipelineError, "'search_paths' is not a ke"
                                                   "y in /Users/ccowart/Develo"
                                                   "pment/mg-scripts/sequence_"
                                                   "processing_pipeline/tests/"
                                                   "data/bad_configuration.jso"
                                                   "n"):
            Pipeline(self.bad_config_file,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        # Pipeline should assert due to an invalid config file path.
        with self.assertRaisesRegex(PipelineError, "does/not/exist/configurati"
                                                   "on.json does not exist."):
            Pipeline(self.invalid_config_file,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        # Pipeline should assert on config_file = None
        with self.assertRaisesRegex(PipelineError, "configuration_file_path ca"
                                                   "nnot be None"):
            Pipeline(None,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        # Pipeline should assert due to invalid_run_id
        with self.assertRaisesRegex(PipelineError, "A run-dir for 'not-sample-"
                                                   "sequence-directory' could "
                                                   "not be found"):
            Pipeline(self.good_config_file,
                     self.invalid_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

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
            # set a range that 211021_A00000_0000_SAMPLE's timestamps must
            # fall within: between 1 and 2 hours.
            cfg['configuration']['pipeline']['younger_than'] = 2
            cfg['configuration']['pipeline']['older_than'] = 1

            # Set directory's timestamp to between one and two hours and
            # verify it is returned by filter_directories_for_time().

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value older than 1 hour ago + 5 min.
            older_than = current_time - (3600 + (5 * 60))
            tp = self.path('211021_A00000_0000_SAMPLE')
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
            utime(self.path('211021_A00000_0000_SAMPLE'), (older_than,
                                                           older_than))

            obs = pipeline.is_within_time_range(tp)
            self.assertEqual(obs, False)

    def test_sample_sheet_validation(self):
        # test successful validation of a good sample-sheet
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_output_file_path, 'my_qiita_id',
                            None)

        msgs, val_sheet = pipeline.validate(self.good_sample_sheet_path)

        msgs = [str(x) for x in msgs]

        # a successful validation should return an empty list of error
        # messages and a sheet object.
        self.assertEqual(msgs, [])
        self.assertIsNotNone(val_sheet)

        # test unsuccessful validation of a bad sample-sheet
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_output_file_path, 'my_qiita_id',
                            None)

        msgs, val_sheet = pipeline.validate(self.bad_sample_sheet_path)

        # an unsuccessful validation should return a list of one or more
        # metapool ErrorMessages and a value of None for val_sheet.
        self.assertIsNotNone(msgs)
        self.assertEqual(str(msgs[0]), 'ErrorMessage: A sample already exists '
                                       'with lane 1 and sample-id EP479894B04')

    def test_generate_sample_information_files(self):
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_output_file_path, 'my_qiita_id',
                            None)

        # assume a validated sample-sheet. Also,
        # generate_sample_information_files calls validate() itself to ensure
        # proper operation.
        paths = partial(pipeline.generate_sample_information_files)
        paths = paths(self.good_sample_sheet_path)

        # confirm files exist in the expected location and with the expected
        # filenames.
        obs = [x.split('sequence_processing_pipeline/')[1] for x in paths]
        exp = ['tests/data/output_dir/NYU_BMS_Melanoma_13059_blanks.tsv',
               'tests/data/output_dir/Feist_11661_blanks.tsv',
               'tests/data/output_dir/Gerwick_6123_blanks.tsv']

        # sort the lists to ensure both are in a fixed order.
        obs.sort()
        exp.sort()

        self.assertEqual(obs, exp)

        # confirm files contain the expected number of lines.
        # This is going to be based on the number of samples named 'BLANK*'
        # in good-sample-sheet.csv.
        exp_lines = {'NYU_BMS_Melanoma_13059_blanks.tsv': 33,
                     'Feist_11661_blanks.tsv': 8,
                     'Gerwick_6123_blanks.tsv': 2}

        exp_first_lines = {
            'NYU_BMS_Melanoma_13059_blanks.tsv': 'BLANK1.1A\t2021-10-21\t193\t'
                                                 'Control\tNegative\tSterile w'
                                                 'ater blank\turban biome\tres'
                                                 'earch facility\tsterile wate'
                                                 'r\tmisc environment\tUSA:CA:'
                                                 'San Diego\tBLANK1.1A\t32.5\t'
                                                 '-117.25\tcontrol blank\tmeta'
                                                 'genome\t256318\tBLANK1.1A\tN'
                                                 'YU_BMS_Melanoma\tTRUE\t'
                                                 'UCSD\tFALSE',
            'Feist_11661_blanks.tsv': 'BLANK.40.12G\t2021-10-21\t193\tControl'
                                      '\tNegative\tSterile water blank\turban '
                                      'biome\tresearch facility\tsterile water'
                                      '\tmisc environment\tUSA:CA:San Diego\tB'
                                      'LANK.40.12G\t32.5\t-117.25\tcontrol bla'
                                      'nk\tmetagenome\t256318\tBLANK.40.12G\t'
                                      'Feist\tTRUE\tUCSD\tFALSE',
            'Gerwick_6123_blanks.tsv': 'BLANK.41.12G\t2021-10-21\t193\tControl'
                                       '\tNegative\tSterile water blank\turban'
                                       ' biome\tresearch facility\tsterile wat'
                                       'er\tmisc environment\tUSA:CA:San Diego'
                                       '\tBLANK.41.12G\t32.5\t-117.25\tcontrol'
                                       ' blank\tmetagenome\t256318\tBLANK.41.1'
                                       '2G\tGerwick\tTRUE\tUCSD\tFALSE'
        }

        exp_last_lines = {
            'NYU_BMS_Melanoma_13059_blanks.tsv': 'BLANK4.4H\t2021-10-21\t193\t'
                                                 'Control\tNegative\tSterile w'
                                                 'ater blank\turban biome\tres'
                                                 'earch facility\tsterile wate'
                                                 'r\tmisc environment\tUSA:CA:'
                                                 'San Diego\tBLANK4.4H\t32.5\t'
                                                 '-117.25\tcontrol blank\tmeta'
                                                 'genome\t256318\tBLANK4.4H\tN'
                                                 'YU_BMS_Melanoma\tTRUE\t'
                                                 'UCSD\tFALSE',
            'Feist_11661_blanks.tsv': 'BLANK.43.12H\t2021-10-21\t193\tControl'
                                      '\tNegative\tSterile water blank\turban'
                                      ' biome\tresearch facility\tsterile wat'
                                      'er\tmisc environment\tUSA:CA:San Diego'
                                      '\tBLANK.43.12H\t32.5\t-117.25\tcontrol'
                                      ' blank\tmetagenome\t256318\tBLANK.43.1'
                                      '2H\tFeist\tTRUE\tUCSD\tFALSE',
            'Gerwick_6123_blanks.tsv': 'BLANK.41.12G\t2021-10-21\t193\tContro'
                                       'l\tNegative\tSterile water blank\turb'
                                       'an biome\tresearch facility\tsterile '
                                       'water\tmisc environment\tUSA:CA:San D'
                                       'iego\tBLANK.41.12G\t32.5\t-117.25\tco'
                                       'ntrol blank\tmetagenome\t256318\tBLAN'
                                       'K.41.12G\tGerwick\tTRUE\tUCSD\t'
                                       'FALSE'
        }

        for some_path in paths:
            some_name = basename(some_path)
            with open(some_path, 'r') as f:
                obs_lines = f.readlines()
                self.assertEqual(len(obs_lines), exp_lines[some_name])
                # confirm that each file contains the expected header.
                header = obs_lines[0].strip()
                self.assertEqual(header, '\t'.join(Pipeline.sif_header))
                # confirm that the first line of each file is as expected.
                obs = obs_lines[1].strip()
                exp = exp_first_lines[some_name]
                self.assertEqual(obs, exp)
                # confirm that the last line of each file is as expected.
                obs = obs_lines[-1].strip()
                exp = exp_last_lines[some_name]
                self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
