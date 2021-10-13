import unittest
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError


class TestJob(unittest.TestCase):
    def test_system_call(self):
        job = Job('/path/to/run_dir', '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], None)

        exp = ('/path/to/run_dir/200nnn_xnnnnn_nnnn_xxxxxxxxxx_1.sh',
               '/path/to/run_dir/200nnn_xnnnnn_nnnn_xxxxxxxxxx_1.out.log',
               '/path/to/run_dir/200nnn_xnnnnn_nnnn_xxxxxxxxxx_1.err.log')
        self.assertEqual(job.generate_job_script_path(), exp)

        self.assertTrue(job._which('ls') in ['/bin/ls', '/usr/bin/ls'])

        with self.assertRaises(PipelineError):
            job._file_check('/does/not/exist')

        self.assertIn(
            './sequence_processing_pipeline/QCCmdGenerator.py',
            job._find_files('./sequence_processing_pipeline'))

        obs = job._system_call('ls sequence_processing_pipeline/'
                               'tests/bin')
        exp = {'stdout': 'bcl-convert\nbcl2fastq\n',
               'stderr': '',
               'return_code': 0}
        self.assertDictEqual(obs, exp)

        with self.assertRaises(PipelineError):
            job._system_call('error')


if __name__ == '__main__':
    unittest.main()
