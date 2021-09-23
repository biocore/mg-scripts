import unittest
from sequence_processing_pipeline.Job import Job
import logging


logging.basicConfig(level=logging.DEBUG)


class TestJob(unittest.TestCase):
    def test_run(self, queue_name, node_count, nprocs, wall_time_limit):
        pass

    def test_directory_check(self, directory_path, create=False):
        pass

    def test_system_call(self, cmd, allow_return_codes=[]):
        job = Job()
        results = job._system_call('ls')

        exp_stdout = {'Job.py', 'FastQC.py', 'PipelineError.py',
                      'dir_locations.txt', '__init__.py',
                      'main.py',
                      'good-bcl-directory', 'TorqueJob.py', 'notes.txt',
                      'SequenceDirectory.py', 'QCJob.py',
                      'metagenomics_pooling_notebook',
                      'BCLConvertJob.py', 'ConvertBCL2FastqJob.py', 'Pipeline.py'}
        exp_stderr = ""
        exp_return_code = 0
        obs_std_out = results['stdout'].split('\n')
        obs_std_out = set([x.strip() for x in obs_std_out if x])
        self.assertEqual(obs_std_out, exp_stdout)
        self.assertEqual(results['stderr'], exp_stderr)
        self.assertEqual(results['return_code'], exp_return_code)

    def test_qsub(self, script_path, qsub_parameters=None, script_parameters=None, wait=True):
        pass
