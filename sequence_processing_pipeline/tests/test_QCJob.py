import shutil
import unittest
from os.path import exists, join, basename, abspath
from functools import partial
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.PipelineError import PipelineError
from os import makedirs
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
import logging

logging.basicConfig(level=logging.INFO)


class TestQCJob(unittest.TestCase):
    def setUp(self):
        # adjustable test_root helps w/testing in different environments.
        test_root = abspath('sequence_processing_pipeline')
        self.path = partial(join, test_root, 'tests', 'data')
        self.sample_sheet_path = self.path('good-sample-sheet.csv')
        self.mmi_db_path = self.path('mmi.db')
        self.sample_run_dir = self.path('MyRunDir')
        self.project_list = ['NYU_BMS_Melanoma_13059', 'Feist_11661',
                             'Gerwick_6123']
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.maxDiff = None

        try:
            shutil.rmtree(self.sample_run_dir)
        except FileNotFoundError:
            # sample_run_dir is created during testing. we don't mind if it
            # doesn't already exist, but any leftovers from previous runs
            # need to be removed.
            pass

        sheet = KLSampleSheet(self.sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        sample_ids = []
        for sample in valid_sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        self.sample_paths = []
        fastq_path = partial(join, self.sample_run_dir, 'Data', 'Fastq')
        for project_name in self.project_list:
            sample_path = fastq_path(project_name)
            makedirs(sample_path, exist_ok=True)
            self.sample_paths.append(sample_path)

            ids = [x[0] for x in
                   filter(lambda c: c[1] == project_name, sample_ids)]

            for id in ids:
                fr_fp = join(sample_path, f'{id}_R1_001.fastq.gz')
                rr_fp = join(sample_path, f'{id}_R2_001.fastq.gz')
                with open(fr_fp, 'w') as f:
                    f.write('This is a forward-read file.')
                with open(rr_fp, 'w') as f:
                    f.write('This is a reverse-read file.')

    def test_split_file_creation(self):
        for sample_path in self.sample_paths:
            qc_job = QCJob(self.sample_run_dir, self.sample_sheet_path,
                           self.mmi_db_path, 'queue_name', 1, 16, 24, '8gb',
                           'fastp', 'minimap2', 'samtools', [],
                           self.qiita_job_id, 30, sample_path)

        # assert that the Torque job files were created and are in the
        # proper location.
        exp = ['QCJob_1.sh', 'QCJob_2.sh', 'QCJob_3.sh']
        exp = set([join(self.sample_run_dir, x) for x in exp])
        obs = set([qc_job.script_paths[proj] for proj in qc_job.script_paths])
        self.assertEqual(obs, exp)

        # compare the expected content for QCJob_1-3.sh with the observed
        # results.
        remove_this = None
        for shell_file in ['QCJob_1.sh', 'QCJob_2.sh', 'QCJob_3.sh']:
            with open(self.path('qcjob_output', shell_file), 'r') as f_exp:
                lines_exp = f_exp.readlines()
                lines_exp = [x.strip() for x in lines_exp]
                with open([x for x in list(exp) if shell_file in x][0],
                          'r') as f_obs:
                    lines_obs = f_obs.readlines()
                    lines_obs = [x.strip() for x in lines_obs]
                    # Technically, remove_this only needs to be calculated
                    # once. Assume there is only one line in the file that
                    # begins w/cd.
                    tmp = [x for x in lines_obs if x.startswith('cd ')][0]
                    # remove the cd, but restore other spaces
                    tmp = ' '.join(tmp.split(' ')[1:])
                    # split off path prepending 'sequence_processing...'
                    tmp = tmp.split('sequence_processing_pipeline/tests')
                    # if tmp2 has '' at element zero, then nothing
                    # prepends 'sequence_processing_pipeline' and nothing
                    # needs to be removed.
                    remove_this = tmp[0] if tmp[0] != '' else None
                    if remove_this:
                        lines_obs = [x.replace(remove_this, '') for x in
                                     lines_obs]

                    self.assertEqual(lines_obs, lines_exp)

        # assert that the array-details files were created and are in the
        # proper location.
        exp = ['split_file_Feist_11661.array-details',
               'split_file_Gerwick_6123.array-details',
               'split_file_NYU_BMS_Melanoma_13059.array-details',
               'split_file_Feist_11661_0', 'split_file_Feist_11661_1',
               'split_file_Feist_11661_2', 'split_file_Feist_11661_3',
               'split_file_Gerwick_6123_0',
               'split_file_NYU_BMS_Melanoma_13059_0',
               'split_file_NYU_BMS_Melanoma_13059_1',
               'split_file_NYU_BMS_Melanoma_13059_2',
               'split_file_NYU_BMS_Melanoma_13059_3']

        exp = set([join(self.sample_run_dir, x) for x in exp])

        for some_path in exp:
            # assert files are in the proper location
            self.assertTrue(exists(some_path) is True)

            # compare the expected content for the paths in exp with the
            # observed results.
            with open(some_path, 'r') as f_obs:
                lines_obs = f_obs.readlines()
                lines_obs = [x.strip() for x in lines_obs]
                # remove_this is saved from above
                if remove_this:
                    lines_obs = [x.replace(remove_this, '') for x in lines_obs]

                with open(self.path('qcjob_output', basename(some_path)),
                          'r') as f_exp:
                    lines_exp = f_exp.readlines()
                    lines_exp = [x.strip() for x in lines_exp]
                    self.assertEqual(lines_obs, lines_exp)

    def test_qcjob_creation(self):
        with self.assertRaises(PipelineError) as e:
            # Pick the first project in the list and generate the path to the
            # Fastq directory.
            fastq_path = partial(join, self.sample_run_dir, 'Data', 'Fastq')
            sample_path = fastq_path(self.project_list[0])

            QCJob(self.sample_run_dir, 'not/path/to/sample/sheet',
                  self.mmi_db_path, 'queue_name', 1, 16, 24, '8gb', 'fastp',
                  'minimap2', 'samtools', [], self.qiita_job_id, 30,
                  sample_path)

        self.assertEqual(str(e.exception), "file 'not/path/to/sample/sheet' "
                                           "does not exist.")

    def test_generate_split_count(self):
        '''
        Extra testing for _generate_split_count() since the number of sample
        projects needed to exercise all of this method would be impractical.
        :return:
        '''
        fastq_path = partial(join, self.sample_run_dir, 'Data', 'Fastq')
        sample_path = fastq_path(self.project_list[0])

        qc_job = QCJob(self.sample_run_dir, self.sample_sheet_path,
                       self.mmi_db_path, 'queue_name', 1, 16, 24, '8gb',
                       'fastp', 'minimap2', 'samtools', [], self.qiita_job_id,
                       30, sample_path)

        input = [-1, 0, 1, 499, 500, 501, 999, 1000, 1001, 1999, 2000, 2001]
        exp_output = [1, 1, 1, 1, 1, 4, 4, 4, 10, 10, 10, 16]
        obs_output = [qc_job._generate_split_count(x) for x in input]

        self.assertEqual(obs_output, exp_output)


if __name__ == '__main__':
    unittest.main()
