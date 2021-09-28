import unittest
from sequence_processing_pipeline.Job import Job
import logging


logging.basicConfig(level=logging.DEBUG)


class TestJob(unittest.TestCase):
    def test_system_call(self):
        job = Job()
        results = job._system_call('ls sequence_processing_pipeline/'
                                   'tests/bin')

        exp_stdout = {'bcl-convert', 'bcl2fastq'}
        exp_stderr = ""
        exp_return_code = 0
        obs_std_out = results['stdout'].split('\n')
        obs_std_out = set([x.strip() for x in obs_std_out if x])
        print(obs_std_out)
        self.assertEqual(obs_std_out, exp_stdout)
        self.assertEqual(results['stderr'], exp_stderr)
        self.assertEqual(results['return_code'], exp_return_code)
