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

        with self.assertRaises(PipelineError) as e:
            job._file_check('/does/not/exist')

        self.assertRegex(str(e.exception),
                         r"^file '/does/not/exist' does not exist.")

        self.assertIn(
            './sequence_processing_pipeline/QCCmdGenerator.py',
            job._find_files('./sequence_processing_pipeline'))

        obs = job._system_call('ls sequence_processing_pipeline/'
                               'tests/bin')
        exp = {'stdout': 'bcl-convert\nbcl2fastq\nfastqc\n',
               'stderr': '',
               'return_code': 0}
        self.assertDictEqual(obs, exp)

        with self.assertRaises(PipelineError) as e:
            job._system_call('error')

        self.assertRegex(str(e.exception), (r"^Execute command-line statement"
                                            r" failure:\nCommand: error\nretur"
                                            r"n code: 127"))


if __name__ == '__main__':
    unittest.main()
