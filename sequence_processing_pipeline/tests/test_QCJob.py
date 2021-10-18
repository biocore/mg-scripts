import shutil
import unittest
from os.path import join
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
        test_root = 'sequence_processing_pipeline'

        self.path = partial(join, test_root, 'tests', 'data')
        self.sample_sheet_path = self.path('good-sample-sheet.csv')
        self.mmi_db_path = self.path('mmi.db')
        self.sample_run_dir = self.path('MyRunDir')

        self.project_list = ['NYU_BMS_Melanoma_13059', 'Feist_11661',
                             'Gerwick_6123']
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'

        self.delete_sample_run_directory()
        self.sample_paths = self.make_sample_run_directory()

    def delete_sample_run_directory(self):
        try:
            shutil.rmtree(self.sample_run_dir)
        except FileNotFoundError:
            # sample_run_dir is created during testing. we don't mind if it
            # doesn't already exist, but any leftovers from previous runs
            # need to be removed.
            pass

    def make_sample_run_directory(self):
        sheet = KLSampleSheet(self.sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        sample_ids = []
        for sample in valid_sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        sample_paths = []
        fastq_path = partial(join, self.sample_run_dir, 'Data', 'Fastq')
        for project_name in self.project_list:
            sample_path = fastq_path(project_name)
            makedirs(sample_path, exist_ok=True)
            sample_paths.append(sample_path)

            ids = [x[0] for x in
                   filter(lambda c: c[1] == project_name, sample_ids)]

            for id in ids:
                fr_fp = join(sample_path, f'{id}_R1_001.fastq.gz')
                rr_fp = join(sample_path, f'{id}_R2_001.fastq.gz')
                with open(fr_fp, 'w') as f:
                    f.write('This is a forward-read file.')
                with open(rr_fp, 'w') as f:
                    f.write('This is a reverse-read file.')

        return sample_paths

    def test_trim_file_creation(self):
        for sample_path in self.sample_paths:
            qc_job = QCJob(self.sample_run_dir, self.sample_sheet_path,
                           self.mmi_db_path, 'queue_name', 1, 16, 24, '8gb',
                           'fastp', 'minimap2', 'samtools', [],
                           self.qiita_job_id, 30, sample_path)

        exp = ['QCJob_1.sh', 'QCJob_2.sh', 'QCJob_3.sh']
        exp = set([join(self.sample_run_dir, x) for x in exp])
        obs = set([qc_job.script_paths[proj] for proj in qc_job.script_paths])

        self.assertEqual(obs, exp)

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


if __name__ == '__main__':
    unittest.main()
