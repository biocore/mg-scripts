import json
import os
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import unittest
from os import makedirs
from os.path import abspath, basename, join
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

        with self.assertRaisesRegex(PipelineError, "required file 'RunInfo.xml"
                                                   "' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)

        # delete RTAComplete.txt and recreate RunInfo.txt file to verify that
        # an Error is raised when only RTAComplete.txt is missing.
        self.delete_rtacomplete_file()
        self.create_runinfo_file()

        with self.assertRaisesRegex(PipelineError, "required file 'RTAComplete"
                                                   ".txt' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)

        # make RunInfo.xml file unreadable and verify that Pipeline object
        # raises the expected Error.
        self.create_rtacomplete_file()
        self.make_runinfo_file_unreadable()

        with self.assertRaisesRegex(PipelineError, "RunInfo.xml is present, bu"
                                                   "t not readable"):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_sample_sheet_path, self.good_output_file_path,
                     self.good_qiita_id, None)
        self.make_runinfo_file_readable()

    def test_creation(self):
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

    def test_sample_sheet_validation(self):
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

        # warning messages may appear before our ErrorMessage but we're not
        # concerned with them at the moment. Extract our ErrorMessage which
        # we'll assume to be the last message.
        exp = str(e.exception).split('\n\n')[-1]
        self.assertEqual(exp, ('ErrorMessage: A sample already exists with lan'
                               'e 1 and sample-id EP479894B04'))

    def test_generate_sample_information_files(self):
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

    def test_get_sample_ids(self):

        exp_sample_ids = ['CDPH-SAL.Salmonella.Typhi.MDL-143',
                          'CDPH-SAL.Salmonella.Typhi.MDL-144',
                          'CDPH-SAL.Salmonella.Typhi.MDL-145',
                          'CDPH-SAL.Salmonella.Typhi.MDL-146',
                          'CDPH-SAL.Salmonella.Typhi.MDL-147',
                          'CDPH-SAL.Salmonella.Typhi.MDL-148',
                          'CDPH-SAL.Salmonella.Typhi.MDL-149',
                          'CDPH-SAL.Salmonella.Typhi.MDL-150',
                          'CDPH-SAL.Salmonella.Typhi.MDL-151',
                          'CDPH-SAL.Salmonella.Typhi.MDL-152',
                          'CDPH-SAL.Salmonella.Typhi.MDL-153',
                          'CDPH-SAL.Salmonella.Typhi.MDL-154',
                          'CDPH-SAL.Salmonella.Typhi.MDL-155',
                          'CDPH-SAL.Salmonella.Typhi.MDL-156',
                          'CDPH-SAL.Salmonella.Typhi.MDL-157',
                          'CDPH-SAL.Salmonella.Typhi.MDL-158',
                          'CDPH-SAL.Salmonella.Typhi.MDL-159',
                          'CDPH-SAL.Salmonella.Typhi.MDL-160',
                          'CDPH-SAL.Salmonella.Typhi.MDL-161',
                          'CDPH-SAL.Salmonella.Typhi.MDL-162',
                          'CDPH-SAL.Salmonella.Typhi.MDL-163',
                          'CDPH-SAL.Salmonella.Typhi.MDL-164',
                          'CDPH-SAL.Salmonella.Typhi.MDL-165',
                          'CDPH-SAL.Salmonella.Typhi.MDL-166',
                          'CDPH-SAL.Salmonella.Typhi.MDL-167',
                          'CDPH-SAL.Salmonella.Typhi.MDL-168',
                          'P21.E.coli.ELI344', 'P21.E.coli.ELI345',
                          'P21.E.coli.ELI347', 'P21.E.coli.ELI348',
                          'P21.E.coli.ELI349', 'P21.E.coli.ELI350',
                          'P21.E.coli.ELI351', 'P21.E.coli.ELI352',
                          'P21.E.coli.ELI353', 'P21.E.coli.ELI354',
                          'P21.E.coli.ELI355', 'P21.E.coli.ELI357',
                          'P21.E.coli.ELI358', 'P21.E.coli.ELI359',
                          'P21.E.coli.ELI361', 'P21.E.coli.ELI362',
                          'P21.E.coli.ELI363', 'P21.E.coli.ELI364',
                          'P21.E.coli.ELI365', 'P21.E.coli.ELI366',
                          'P21.E.coli.ELI367', 'P21.E.coli.ELI368',
                          'P21.E.coli.ELI369', 'stALE.E.coli.A1.F21.I1.R1',
                          'stALE.E.coli.A2.F21.I1.R1',
                          'stALE.E.coli.A3.F18.I1.R1',
                          'stALE.E.coli.A3.F40.I1.R1',
                          'stALE.E.coli.A4.F21.I1.R1',
                          'stALE.E.coli.A4.F21.I1.R2',
                          'stALE.E.coli.A4.F42.I1.R1',
                          'stALE.E.coli.A5.F21.I1.R1',
                          'stALE.E.coli.A5.F42.I1.R1',
                          'stALE.E.coli.A6.F21.I1.R1',
                          'stALE.E.coli.A6.F43.I1.R1',
                          'stALE.E.coli.A7.F21.I1.R1',
                          'stALE.E.coli.A7.F42.I1.R1',
                          'stALE.E.coli.A8.F20.I1.R1',
                          'stALE.E.coli.A8.F42.I1.R1',
                          'stALE.E.coli.A9.F21.I1.R1',
                          'stALE.E.coli.A9.F44.I1.R1',
                          'stALE.E.coli.A10.F21.I1.R1',
                          'stALE.E.coli.A10.F43.I1.R1',
                          'stALE.E.coli.A10.F131.I1.R1',
                          'stALE.E.coli.A11.F21.I1.R1',
                          'stALE.E.coli.A11.F43.I1.R1',
                          'stALE.E.coli.A11.F119.I1.R1',
                          'stALE.E.coli.A12.F21.I1.R1',
                          'stALE.E.coli.A12.F43.I1.R1',
                          'stALE.E.coli.A12.F136.I1.R1',
                          'stALE.E.coli.A13.F20.I1.R1',
                          'stALE.E.coli.A13.F42.I1.R1',
                          'stALE.E.coli.A13.F121.I1.R1',
                          'stALE.E.coli.A14.F20.I1.R1',
                          'stALE.E.coli.A14.F42.I1.R1',
                          'stALE.E.coli.A14.F133.I1.R1',
                          'stALE.E.coli.A15.F21.I1.R1',
                          'stALE.E.coli.A15.F42.I1.R1',
                          'stALE.E.coli.A15.F117.I1.R1',
                          'stALE.E.coli.A16.F20.I1.R1',
                          'stALE.E.coli.A16.F42.I1.R1',
                          'stALE.E.coli.A16.F134.I1.R1',
                          'stALE.E.coli.A17.F21.I1.R1',
                          'stALE.E.coli.A17.F118.I1.R1',
                          'stALE.E.coli.A18.F18.I1.R1',
                          'stALE.E.coli.A18.F39.I1.R1',
                          'stALE.E.coli.A18.F130.I1.R1', '3A', '4A',
                          'BLANK.40.12G', 'BLANK.40.12H',
                          'Pputida.JBEI.HGL.Pputida.107.BP6',
                          'Pputida.JBEI.HGL.Pputida.108.BP7',
                          'Pputida.JBEI.HGL.Pputida.109.BP8',
                          'Pputida.JBEI.HGL.Pputida.110.M2',
                          'Pputida.JBEI.HGL.Pputida.111.M5',
                          'Pputida.TALE.HGL.Pputida.112',
                          'Pputida.TALE.HGL.Pputida.113',
                          'Pputida.TALE.HGL.Pputida.114',
                          'Pputida.TALE.HGL.Pputida.115',
                          'Pputida.TALE.HGL.Pputida.116',
                          'Pputida.TALE.HGL.Pputida.117',
                          'Pputida.TALE.HGL.Pputida.118',
                          'Pputida.TALE.HGL.Pputida.119',
                          'Pputida.TALE.HGL.Pputida.120',
                          'Pputida.TALE.HGL.Pputida.121',
                          'Pputida.TALE.HGL.Pputida.122',
                          'Pputida.TALE.HGL.Pputida.123',
                          'Pputida.TALE.HGL.Pputida.124',
                          'Pputida.TALE.HGL.Pputida.125',
                          'Pputida.TALE.HGL.Pputida.126',
                          'Pputida.TALE.HGL.Pputida.127',
                          'Pputida.TALE.HGL.Pputida.128',
                          'Pputida.TALE.HGL.Pputida.129',
                          'Pputida.TALE.HGL.Pputida.130',
                          'Pputida.TALE.HGL.Pputida.131',
                          'Pputida.TALE.HGL.Pputida.132',
                          'Pputida.TALE.HGL.Pputida.133',
                          'Pputida.TALE.HGL.Pputida.134',
                          'Pputida.TALE.HGL.Pputida.135',
                          'Pputida.TALE.HGL.Pputida.136',
                          'Pputida.TALE.HGL.Pputida.137',
                          'Pputida.TALE.HGL.Pputida.138',
                          'Pputida.TALE.HGL.Pputida.139',
                          'Pputida.TALE.HGL.Pputida.140',
                          'Pputida.TALE.HGL.Pputida.141',
                          'Pputida.TALE.HGL.Pputida.142',
                          'Pputida.TALE.HGL.Pputida.143',
                          'Pputida.TALE.HGL.Pputida.144',
                          'Pputida.PALE.HGL.Pputida.145',
                          'Pputida.PALE.HGL.Pputida.146',
                          'Pputida.PALE.HGL.Pputida.147',
                          'Pputida.PALE.HGL.Pputida.148',
                          'Pputida.PALE.HGL.Pputida.149',
                          'Pputida.PALE.HGL.Pputida.150',
                          'Pputida.PALE.HGL.Pputida.151',
                          'Pputida.PALE.HGL.Pputida.152',
                          'Pputida.PALE.HGL.Pputida.153',
                          'Pputida.PALE.HGL.Pputida.154',
                          'Pputida.PALE.HGL.Pputida.155',
                          'Pputida.PALE.HGL.Pputida.156',
                          'Pputida.PALE.HGL.Pputida.157',
                          'Pputida.PALE.HGL.Pputida.158',
                          'Pputida.PALE.HGL.Pputida.159',
                          'Pputida.PALE.HGL.Pputida.160',
                          'Pputida.PALE.HGL.Pputida.161',
                          'Pputida.PALE.HGL.Pputida.162',
                          'Pputida.PALE.HGL.Pputida.163',
                          'Pputida.PALE.HGL.Pputida.164',
                          'Pputida.PALE.HGL.Pputida.165',
                          'Pputida.PALE.HGL.Pputida.166',
                          'Pputida.PALE.HGL.Pputida.167',
                          'Pputida.PALE.HGL.Pputida.168',
                          'Pputida.PALE.HGL.Pputida.169',
                          'Pputida.PALE.HGL.Pputida.170',
                          'Pputida.PALE.HGL.Pputida.171',
                          'Pputida.PALE.HGL.Pputida.172',
                          'Pputida.PALE.HGL.Pputida.173',
                          'Pputida.PALE.HGL.Pputida.174',
                          'Pputida.PALE.HGL.Pputida.175',
                          'Pputida.PALE.HGL.Pputida.176',
                          'JM-Metabolic.GN0.2005', 'JM-Metabolic.GN0.2007',
                          'JM-Metabolic.GN0.2009', 'JM-Metabolic.GN0.2094',
                          'JM-Metabolic.GN0.2099', 'JM-Metabolic.GN0.2148',
                          'JM-Metabolic.GN0.2165', 'JM-Metabolic.GN0.2169',
                          'JM-Metabolic.GN0.2172', 'JM-Metabolic.GN0.2175',
                          'JM-Metabolic.GN0.2183', 'JM-Metabolic.GN0.2215',
                          'JM-Metabolic.GN0.2254', 'JM-Metabolic.GN0.2277',
                          'JM-Metabolic.GN0.2290', 'JM-Metabolic.GN0.2337',
                          'JM-Metabolic.GN0.2317', 'JM-Metabolic.GN0.2354',
                          'JM-Metabolic.GN0.2375', 'JM-Metabolic.GN0.2380',
                          'JM-Metabolic.GN0.2393', 'JM-Metabolic.GN0.2404',
                          '5B', '6A', 'BLANK.41.12G', 'BLANK.41.12H',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.4.14',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.4.23',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.4.48',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.6.21',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.6.35',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.10.13',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.10.28',
                          'Deoxyribose.PALE.ALE.MG1655.BOP27.10.51',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.18.19',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.18.59',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.18.35',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.20.16',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.20.43',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.20.71',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.22.16',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.22.28',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.22.52',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.24.9',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.24.24',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.24.52',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.26.6',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.26.27',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.26.69',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.28.13',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.28.28',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.28.53',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.30.7',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.30.22',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.30.60',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.32.6',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.32.20',
                          'Deoxyribose.PALE.ALE.MG1655.Lib4.32.56',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.24',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.57',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.1.69',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.23',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.50',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.3.61',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.22',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.36',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.5.46',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.23',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.41',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.7.51',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.25',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.58',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.17.64',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.25',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.55',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.19.63',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.23',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.46',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.21.51',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.25',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.49',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.29.57',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.24',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.42',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.31.62',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.21',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.41',
                          'AB5075.AZM.TALE.in.MHB.A.baumannii.AB5075.WT.33.50',
                          'JM-Metabolic.GN02514', 'JM-Metabolic.GN02529',
                          'JM-Metabolic.GN02531', 'JM-Metabolic.GN02567',
                          'JM-Metabolic.GN02590', 'JM-Metabolic.GN02657',
                          'JM-Metabolic.GN02748', 'JM-Metabolic.GN02766',
                          'JM-Metabolic.GN02769', 'JM-Metabolic.GN02787',
                          'JM-Metabolic.GN03132', 'JM-Metabolic.GN03218',
                          'JM-Metabolic.GN03252', 'JM-Metabolic.GN03409',
                          'JM-Metabolic.GN04014', 'JM-Metabolic.GN04094',
                          'JM-Metabolic.GN04255', 'JM-Metabolic.GN04306',
                          'JM-Metabolic.GN04428', 'JM-Metabolic.GN04488',
                          'JM-Metabolic.GN04540', 'JM-Metabolic.GN04563',
                          'JM-Metabolic.GN04612', 'JM-Metabolic.GN04665',
                          'JM-Metabolic.GN04682', 'JM-Metabolic.GN05002',
                          'JM-Metabolic.GN05109', 'JM-Metabolic.GN05128',
                          'JM-Metabolic.GN05367', 'JM-Metabolic.GN05377',
                          '7A', '8A', 'BLANK.42.12G', 'BLANK.42.12H',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0326',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0327',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0328',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0329',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0330',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0352',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0353',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0354',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0355',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0356',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0357',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0364',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0366',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0367',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0368',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0369',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0370',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0371',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0372',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0373',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0374',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0375',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0376',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0377',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0378',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0380',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0381',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0382',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0383',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0384',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0385',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0386',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0387',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0388',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0389',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0390',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0391',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0392',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0393',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0394',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0395',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0396',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0397',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0398',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0399',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0400',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0401',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0402',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0403',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0404',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0405',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0406',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0407',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0408',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0409',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0417',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0418',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0419',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0420',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0421',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0473',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0474',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0483',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0484',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0485',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0486',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0516',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0517',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0518',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0519',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0520',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0521',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0522',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0523',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0524',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-B0525',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R08624',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R08704',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R10727',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11044',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11078',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11101',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11102',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11103',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11135',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11153',
                          'JM-MEC.Staphylococcus.aureusstrain.BERTI-R11154',
                          'JM-Metabolic.GN02424', 'JM-Metabolic.GN02446',
                          'JM-Metabolic.GN02449', 'JM-Metabolic.GN02487',
                          'JM-Metabolic.GN02501', 'ISB', 'GFR',
                          'BLANK.43.12G',
                          'BLANK.43.12H', 'RMA.KHP.rpoS.Mage.Q97D',
                          'RMA.KHP.rpoS.Mage.Q97L', 'RMA.KHP.rpoS.Mage.Q97N',
                          'RMA.KHP.rpoS.Mage.Q97E', 'JBI.KHP.HGL.021',
                          'JBI.KHP.HGL.022', 'JBI.KHP.HGL.023',
                          'JBI.KHP.HGL.024', 'JBI.KHP.HGL.025',
                          'JBI.KHP.HGL.026', 'JBI.KHP.HGL.027',
                          'JBI.KHP.HGL.028.Amitesh.soxR',
                          'JBI.KHP.HGL.029.Amitesh.oxyR',
                          'JBI.KHP.HGL.030.Amitesh.soxR.oxyR',
                          'JBI.KHP.HGL.031.Amitesh.rpoS', 'BLANK1.1A',
                          'BLANK1.1B', 'BLANK1.1C', 'BLANK1.1D', 'BLANK1.1E',
                          'BLANK1.1F', 'BLANK1.1G', 'BLANK1.1H', 'AP581451B02',
                          'EP256645B01', 'EP112567B02', 'EP337425B01',
                          'LP127890A01', 'EP159692B04', 'EP987683A01',
                          'AP959450A03', 'SP464350A04', 'C9', 'ep256643b01',
                          'EP121011B01', 'AP616837B04', 'SP506933A04',
                          'EP159695B01', 'EP256644B01', 'SP511289A02',
                          'EP305735B04', 'SP415030A01', 'AP549681B02',
                          'AP549678B01', 'EP260544B04', 'EP202452B01',
                          'EP282276B04', 'SP531696A04', 'SP515443A04',
                          'SP515763A04', 'EP184255B04', 'SP503615A02',
                          'EP260543B04', 'EP768748A04', 'AP309872B03',
                          'AP568785B04', 'EP721390A04', 'EP940013A01',
                          'EP291979B04', 'EP182065B04', 'EP128904B02',
                          'EP915769A04', 'SP464352A03', 'SP365864A04',
                          'SP511294A04', 'EP061002B01', 'SP410793A01',
                          'SP232077A04', 'EP128910B01', 'AP531397B04',
                          'EP043583B01', 'EP230245B01', 'EP606652B04',
                          'EP207041B01', 'EP727972A04', 'EP291980B04',
                          'EP087938B02', 'SP471496A04', 'SP573823A04',
                          'EP393718B01', 'SP612496A01', 'EP032410B02',
                          'EP073216B01', 'EP410046B01', 'SP561451A04',
                          'EP320438B01', 'SP612495A04', 'EP446604B03',
                          'EP446602B01', 'EP182243B02', 'EP333541B04',
                          'EP238034B01', 'AP298002B02', 'EP455759B04',
                          'EP207042B04', 'LP128479A01', 'LP128476A01',
                          'EP316863B03', 'C20', 'lp127896a01', 'SP491907A02',
                          'EP182060B03', 'EP422407B01', 'SP573859A04',
                          'SP584547A02', 'EP182346B04', 'AP668631B04',
                          'EP451428B04', 'LP128538A01', 'SP490298A02',
                          'SP573860A01', 'EP032412B02', 'EP163771B01',
                          'LP169879A01', 'EP729433A02', 'EP447940B04',
                          'SP584551A08', 'EP216516B04', 'EP023808B02',
                          'BLANK2.2A', 'BLANK2.2B', 'BLANK2.2C', 'BLANK2.2D',
                          'BLANK2.2E', 'BLANK2.2F', 'BLANK2.2G', 'BLANK2.2H',
                          'SP573843A04', 'EP683835A01', 'SP573824A04',
                          'SP335002A04', 'SP478193A02', 'SP232311A04',
                          'SP415021A02', 'SP231630A02', 'SP641029A02',
                          'SP232310A04', 'EP617442B01', 'EP587478B04',
                          'EP447928B04', 'EP587475B04', 'EP675042B01',
                          'EP554513B02', 'EP702221B04', 'AP568787B02',
                          'EP054632B01', 'EP121013B01', 'EP649418A02',
                          'EP573313B01', 'LP154981A01', 'AP470859B01',
                          'LP154986A01', 'AP732307B04', 'EP533426B03',
                          'EP587476B04', 'AP696363B02', 'EP587477B04',
                          'SP683466A02', 'EP554518B04', 'EP533429B04',
                          'EP431570B01', 'EP202095B04', 'EP504030B04',
                          'EP207036B01', 'EP393717B01', 'SP491898A02',
                          'EP484973B04', 'EP479794B02', 'EP554515B04',
                          'SP631994A04', 'EP921593A04', 'AP787247B04',
                          'EP090129B04', 'EP447975B02', 'EP212214B01',
                          'EP410042B01', 'SP404409A02', 'SP247340A04',
                          'AP029018B01', 'EP872341A01', 'AP062219B03',
                          'EP790020A02', 'EP808112A04', 'SP404403A02',
                          'EP073160B01', 'EP012991B03', 'SP317297A02',
                          'EP656055A04', 'EP649623A01', 'EP790019A01',
                          'SP257519A04', 'EP808104A01', 'EP808106A01',
                          'SP231629A02', 'EP675044A01', 'EP657260A01',
                          'EP808110A04', 'AP032413B04', 'EP843906A04',
                          'AP173305B04', 'SP231628A02', 'AP173301B04',
                          'SP404405A02', 'EP649653A04', 'EP718687A04',
                          'AP905750A02', 'EP738468A01', 'C6', 'EP890157A02',
                          'SP353893A02', 'EP944059A02', 'EP970005A01',
                          'EP927461A04', 'EP808111A03', 'EP927459A04',
                          'SP317293A02', 'SP235186A04', 'SP399724A04',
                          'EP738469A01', 'SP284095A03', 'C5', 'EP337325B04',
                          'EP759450A04', 'BLANK3.3A', 'BLANK3.3B', 'BLANK3.3C',
                          'BLANK3.3D', 'BLANK3.3E', 'BLANK3.3F', 'BLANK3.3G',
                          'BLANK3.3H', 'AP006367B02', 'EP929277A02',
                          'AP324642B04', 'EP786631A04', 'EP657385A04',
                          'SP235189A01', 'EP448041B04', 'SP231631A02',
                          'SP280481A02', 'AP032412B04', 'EP649737A03',
                          'AP967057A04', 'EP876243A04', 'SP229387A04',
                          'EP667743A04', 'SP246941A01', 'AP745799A04',
                          'SP205732A02', 'SP230382A04', 'SP230380A02',
                          'SP230381A01', 'SP205754A01', 'EP606662B04',
                          'AP780167B02', 'EP447927B04', 'C18', 'LP191039A01',
                          'EP606663B04', 'EP573296B01', 'EP447926B04',
                          'LP127767A01', 'EP479266B04', 'LP128543A01',
                          'EP479270B03', 'EP921594A04', 'EP554501B04',
                          'EP542577B04', 'EP487995B04', 'EP542578B04',
                          'EP573310B01', 'EP244366B01', 'EP533389B03',
                          'EP244360B01', 'AP911328B01', 'AP481403B02',
                          '22.001.801.552.503.00', 'EP372981B04',
                          'EP447929B04', 'SP573849A04', 'SP577399A02',
                          'EP606656B03', 'LP166715A01', 'AP668628B04',
                          'C14', 'EP446610B02', 'EP339061B02', 'SP681591A04',
                          'EP393712B02', 'EP410041B01', 'SP453872A01',
                          '22.001.710.503.791.00',
                          'LP128540A01', 'EP339053B02', 'EP617443B01',
                          'EP190307B01', 'AP795068B04', 'LP128541A01',
                          'EP584756B04', 'SP284096A02', 'EP431562B04',
                          'EP685640B01', 'EP339059B02', 'EP431575B01',
                          'EP379938B01', 'EP529635B02', 'EP554506B04',
                          'EP455757B04', 'SP491900A02', 'LP196272A01',
                          'SP704319A04', 'EP617441B01', 'AP687591B04',
                          'SP640978A02', 'EP981129A02', 'EP455763B04',
                          'EP339057B02', 'SP491897A02', 'EP980752B04',
                          'LP128539A01', 'EP996831B04', 'EP273332B04',
                          'EP483291B04', 'EP393715B01', 'EP617440B01',
                          'EP729434A01', 'SP645141A03', 'BLANK4.4A',
                          'BLANK4.4B', 'BLANK4.4C', 'BLANK4.4D', 'BLANK4.4E',
                          'BLANK4.4F', 'BLANK4.4G', 'BLANK4.4H', 'SP232114A04',
                          'EP393714B01', 'EP533388B01', 'EP724905B01',
                          'EP282108B01', 'EP282107B01', 'EP001625B01',
                          'EP073209B02', 'SP232079A01', 'EP772145A02',
                          'AP771472A04', 'AP223470B01', 'SP404412A02',
                          'EP772143A02', 'SP408629A01', 'EP749735A07',
                          'EP846485A01', 'EP808109A01', 'SP416130A04',
                          'EP882752A01', 'AP953594A02', 'AP046324B02',
                          'AP891020A04', 'EP790023A01', 'EP657386A01',
                          'EP805337A01', 'EP927458A04', 'AP173299B04',
                          'EP768164A02', 'EP886422A01', 'AP103463B01',
                          'AP744361A02', 'AP065292B01', 'SP257517A04',
                          'EP790021A04', 'EP675075A04', 'SP388683A02',
                          'SP232309A01', 'EP899038A04', 'EP636802A01',
                          'AP046327B02', 'EP905975A04', 'SP410796A02',
                          'EP784608A01', 'EP808105A01', 'SP331134A04',
                          'EP718688A01', 'SP232270A02', 'EP970001A01',
                          'EP001624B01', 'EP868682A01', 'EP927462A02', 'C3',
                          'EP890158A02', 'EP023801B04', 'EP400447B04',
                          'EP385379B01', 'EP385387B01', 'EP385384B01',
                          'SP754514A04', 'SP415025A01', 'SP415023A02',
                          'EP400448B04', 'EP479894B04']
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.good_output_file_path, self.good_qiita_id,
                            None)

        obs = pipeline.get_sample_ids()
        self.assertEqual(sorted(obs), sorted(exp_sample_ids))

    def test_get_project_info(self):
        exp_proj_info = [
            {'project_name': 'NYU_BMS_Melanoma_13059', 'qiita_id': '13059'},
            {'project_name': 'Feist_11661', 'qiita_id': '11661'},
            {'project_name': 'Gerwick_6123', 'qiita_id': '6123'}]

        exp_project_names = ['NYU_BMS_Melanoma_13059', 'Feist_11661',
                             'Gerwick_6123']

        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_sample_sheet_path,
                            self.good_output_file_path, self.good_qiita_id,
                            None)

        obs_proj_info = pipeline.get_project_info()
        obs_project_names = []
        for d in obs_proj_info:
            obs_project_names.append(d['project_name'])

        self.assertEqual(sorted(obs_project_names), sorted(exp_project_names))

        for exp_d in exp_proj_info:
            for obs_d in obs_proj_info:
                if obs_d['project_name'] == exp_d['project_name']:
                    self.assertDictEqual(obs_d, exp_d)
                    break


if __name__ == '__main__':
    unittest.main()
