from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from os.path import join
from os import makedirs
import unittest
import shutil
import re


class TestGenPrepFileJob(unittest.TestCase):
    def setUp(self):
        self.package_root = 'sequence_processing_pipeline'
        self.qiita_job_id = 'b197f317-1c06-4619-9af3-65721149c1e8'
        self.working_directory_root = join(self.package_root,
                                           self.qiita_job_id)
        try:
            shutil.rmtree(self.working_directory_root)
        except FileNotFoundError:
            # Clean up test directory just in case
            pass
        makedirs(self.working_directory_root)
        self.run_id = '210518_A00953_0305_TEST'
        self.run_dir = join(self.working_directory_root, self.run_id)
        self.convert_job_path = join(self.run_dir, 'ConvertJob')
        makedirs(join(self.convert_job_path, 'Reports'), exist_ok=True)
        self.qc_job_path = join(self.run_dir, 'QCJob')
        self.project_list = ['Project1']
        makedirs(join(self.qc_job_path, self.project_list[0],
                      'filtered_sequences'))
        makedirs(join(self.qc_job_path, self.project_list[0],
                      'fastp_reports_dir', 'json'))

    def tearDown(self):
        shutil.rmtree(self.working_directory_root)

    def test_creation(self):
        self.maxDiff = None
        sample_sheet_path = join(self.package_root, 'tests', 'data',
                                 'good-sample-sheet.csv')

        job = GenPrepFileJob(self.run_dir,
                             self.convert_job_path,
                             self.qc_job_path,
                             join(self.working_directory_root, 'OutputPath'),
                             sample_sheet_path,
                             'seqpro',
                             self.project_list,
                             [],
                             'abcdabcdabcdabcdabcdabcdabcdabcd')
        results = job._system_call(f'find {self.run_dir}')
        lines = results['stdout'].split('\n')
        lines = [re.sub(r'^.*?\/sequence_processing_pipeline\/', '', x)
                 for x in lines]
        obs = [x for x in lines if x]
        obs.sort()

        exp = ['b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'ConvertJob',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'ConvertJob/Reports',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1/fastp_reports_dir',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1/fastp_reports_dir/json',
               'b197f317-1c06-4619-9af3-65721149c1e8/210518_A00953_0305_TEST/'
               'QCJob/Project1/filtered_sequences']

        self.assertEquals(obs, exp)


if __name__ == '__main__':
    unittest.main()
