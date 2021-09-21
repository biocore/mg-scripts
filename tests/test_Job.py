import unittest
from sequence_processing_pipeline.Job import lsJob


class TestJob(unittest.TestCase):
    def test_job(self):
        # Job is a minimal base-class. Use lsJob to test basic Job
        # functionality.
        job = lsJob('sequence_processing_pipeline')
        results = job.run()
        exp_stdout = {'Job.py', 'FastQC.py', 'PipelineError.py',
                      'dir_locations.txt', '__init__.py',
                      'main.py',
                      'good-bcl-directory', 'TorqueJob.py', 'notes.txt',
                      'SequenceDirectory.py', 'HumanFilterJob.py',
                      'metagenomics_pooling_notebook',
                      'BCLConvertJob.py', 'BCL2FASTQJob.py', 'Pipeline.py'}
        exp_stderr = ""
        exp_return_code = 0
        obs_std_out = results['stdout'].split('\n')
        obs_std_out = set([x.strip() for x in obs_std_out if x])
        self.assertEqual(obs_std_out, exp_stdout)
        self.assertEqual(results['stderr'], exp_stderr)
        self.assertEqual(results['return_code'], exp_return_code)



