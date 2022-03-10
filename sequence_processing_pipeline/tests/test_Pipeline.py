import json
import os
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import unittest
from os import utime, makedirs
from os.path import abspath, basename, join
from copy import deepcopy
from time import time
from functools import partial
import re


class TestPipeline(unittest.TestCase):
    def setUp(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.bad_config_file = self.path('bad_configuration.json')
        self.invalid_config_file = 'does/not/exist/configuration.json'
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        # good_qiita_id is randomly-generated and does not match any known
        # existing qiita job_id.
        self.good_qiita_id = '077c4da8-74eb-4184-8860-0207f53623be'
        self.invalid_run_id = 'not-sample-sequence-directory'
        self.good_output_file_path = self.path('output_dir')
        makedirs(self.good_output_file_path, exist_ok=True)
        self.maxDiff = None
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.bad_sample_sheet_path = self.path('duplicate_sample-sample-sheet'
                                               '.csv')
        self.good_run_dir = self.path(self.good_run_id)
        self.runinfo_file = self.path(self.good_run_id, 'RunInfo.xml')
        self.rtacomplete_file = self.path(self.good_run_id, 'RTAComplete.txt')

        # most of the tests here were written with the assumption that these
        # files already exist.
        self.create_runinfo_file()
        self.create_rtacomplete_file()

        # read good configuration file at initialization to avoid putting
        # most of the test code within 'with open()' expressions.
        with open(self.good_config_file, 'r') as f:
            self.good_config = json.load(f)

        # create a fake age for self.good_run_dir by taking the current time
        # and subtracting from it a value in seconds that's in the middle of
        # the allowable age-range.
        config = self.good_config['configuration']['pipeline']
        valid_age = ((config['younger_than'] - config[
                      'older_than']) / 2) * 3600

        self.new_timestamp = time() - valid_age

    def reset_run_directory_timestamp(self, timestamp=None):
        '''
        Reset the timestamp of self.good_run_dir.
        :param timestamp: Reset to supplied value or known good value if None.
        :return:
        '''
        ts = timestamp if timestamp is not None else self.new_timestamp
        utime(self.good_run_dir, (ts, ts))

    def tearDown(self):
        # Pipeline is now the only class aware of these files, hence they
        # can be deleted at the end of testing.
        self.delete_runinfo_file()
        self.delete_rtacomplete_file()

    def make_runinfo_file_unreadable(self):
        os.chmod(self.runinfo_file, 0o000)

    def make_runinfo_file_readable(self):
        os.chmod(self.runinfo_file, 0o777)

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

        # reset the timestamp on self.good_run_dir because it was just
        # created.
        self.reset_run_directory_timestamp()

        with self.assertRaisesRegex(PipelineError, "required file 'RunInfo.xml"
                                                   "' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)

        # delete RTAComplete.txt and recreate RunInfo.txt file to verify that
        # an Error is raised when only RTAComplete.txt is missing.
        self.delete_rtacomplete_file()
        self.create_runinfo_file()

        # with the above changes we need to adjust the timestamp again.
        self.reset_run_directory_timestamp()

        with self.assertRaisesRegex(PipelineError, "required file 'RTAComplete"
                                                   ".txt' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)

        # make RunInfo.xml file unreadable and verify that Pipeline object
        # raises the expected Error.
        self.create_rtacomplete_file()
        self.make_runinfo_file_unreadable()

        # with the above changes we need to adjust the timestamp again.
        self.reset_run_directory_timestamp()

        with self.assertRaisesRegex(PipelineError, "RunInfo.xml is present, bu"
                                                   "t not readable"):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)
        self.make_runinfo_file_readable()

    def test_creation(self):
        # reset the timestamp on self.good_run_dir because it was just
        # created.
        self.reset_run_directory_timestamp()

        # Pipeline should assert due to config_file
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.bad_config_file,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.good_output_file_path,
                     self.good_qiita_id,
                     None)

        msg = re.sub(r'not a key in .*?/sequence_processing_pipeline',
                     r'not a key in sequence_processing_pipeline',
                     str(e.exception))
        self.assertEqual(msg, "'search_paths' is not a key in "
                              "sequence_processing_pipeline/tests"
                              "/data/bad_configuration.json")

        # Pipeline should assert due to an invalid config file path.
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.invalid_config_file,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.good_output_file_path,
                     self.good_qiita_id,
                     None)

        self.assertEqual(str(e.exception), 'does/not/exist/configuration.json '
                                           'does not exist.')

        # Pipeline should assert on config_file = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(None,
                     self.good_run_id,
                     self.good_sample_sheet_path,
                     self.good_output_file_path,
                     self.good_qiita_id,
                     None)

        self.assertEqual(str(e.exception), 'configuration_file_path cannot be '
                                           'None')

        # Pipeline should assert due to invalid_run_id
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     self.invalid_run_id,
                     self.good_sample_sheet_path,
                     self.good_output_file_path,
                     self.good_qiita_id,
                     None)

        self.assertEqual(str(e.exception), "A run-dir for 'not-sample-sequence"
                                           "-directory' could not be found")

        # Pipeline should assert on run_id = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     None,
                     self.good_sample_sheet_path,
                     self.good_output_file_path,
                     self.good_qiita_id,
                     None)

        with open(join('sequence_processing_pipeline', 'configuration.json'),
                  'r') as f:
            cfg = json.load(f)
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = -1
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_sample_sheet_path,
                         self.good_output_file_path,
                         self.good_qiita_id,
                         cpy_cfg)

            self.assertEqual(str(e.exception), 'older_than and younger_than '
                                               'cannot be less than zero.')
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['older_than'] = -1
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_sample_sheet_path,
                         self.good_output_file_path,
                         self.good_qiita_id,
                         cpy_cfg)
            self.assertEqual(str(e.exception), 'older_than and younger_than '
                                               'cannot be less than zero.')
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = 20
                cpy_cfg['configuration']['pipeline']['older_than'] = 30
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_sample_sheet_path,
                         self.good_output_file_path,
                         self.good_qiita_id,
                         cpy_cfg)
            self.assertEqual(str(e.exception), 'older_than cannot be equal to '
                                               'or less than younger_than.')

    def test_filter_directories_for_time(self):
        cfg = deepcopy(self.good_config)

        # set a range that 211021_A00000_0000_SAMPLE's timestamps must
        # fall within: between 1 and 2 hours.
        cfg['configuration']['pipeline']['younger_than'] = 2
        cfg['configuration']['pipeline']['older_than'] = 1

        # Set directory's timestamp to between one and two hours and
        # verify it is returned by filter_directories_for_time().

        # get the current time in seconds since the epoch.
        current_time = time()
        # create an epoch time value older than 1 hour ago + 5 min.
        older_than = current_time - 3900
        self.reset_run_directory_timestamp(timestamp=older_than)

        try:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, cfg)
        except PipelineError as e:
            self.fail(("test_filter_directories_for_time failed w/PipelineEr"
                       f"ror: {e.message}"))

        # Set directory's timestamp to just older than two hours and verify it
        # is not returned by filter_directories_for_time().

        # get the current time in seconds since the epoch.
        current_time = time()
        # create an epoch time value older than 2 hour ago + 5 min.
        older_than = current_time - 7500
        self.reset_run_directory_timestamp(timestamp=older_than)

        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, cfg)
        self.assertEqual(str(e.exception), ('sequence_processing_pipeline/tes'
                                            'ts/data/211021_A00000_0000_SAMPL'
                                            'E is too old.'))

        # Set directory's timestamp to just under one hour and
        # verify it is not returned by filter_directories_for_time().

        # get the current time in seconds since the epoch.
        current_time = time()
        # create an epoch time value younger than 1 hour ago.
        older_than = current_time - 3300
        self.reset_run_directory_timestamp(timestamp=older_than)

        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, cfg)
        self.assertEqual(str(e.exception), ('sequence_processing_pipeline/tes'
                                            'ts/data/211021_A00000_0000_SAMPL'
                                            'E is too young.'))

    def test_sample_sheet_validation(self):
        # reset the timestamp on self.good_run_dir because it was just
        # created.
        self.reset_run_directory_timestamp()

        # test successful validation of a good sample-sheet.
        # if self.good_sample_sheet_path points to a bad sample-sheet, then
        # Pipeline would raise a PipelineError w/warnings and error messages
        # contained w/in its 'message' member.
        try:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)
        except PipelineError as e:
            self.fail(("test_filter_directories_for_time failed w/PipelineEr"
                       f"ror: {e.message}"))

        # test unsuccessful validation of a bad sample-sheet
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file, self.good_run_id,
                     self.bad_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)
        self.assertEqual(str(e.exception), ('ErrorMessage: A sample already e'
                                            'xists with lane 1 and sample-id '
                                            'EP479894B04'))

    def test_generate_sample_information_files(self):
        # reset the timestamp on self.good_run_dir because it was just
        # created.
        self.reset_run_directory_timestamp()

        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.good_output_file_path, self.good_qiita_id,
                            None)

        paths = pipeline.generate_sample_information_files()

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
