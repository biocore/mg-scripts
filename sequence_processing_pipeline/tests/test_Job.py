import unittest
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import abspath, join
from functools import partial


class TestJob(unittest.TestCase):
    def test_system_call(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

        job = Job(self.path('sample-sequence-directory'),
                  self.path('my_output_dir'), '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], None)

        self.assertTrue(job._which('ls') in ['/bin/ls', '/usr/bin/ls'])

        with self.assertRaises(PipelineError) as e:
            job._file_check('/does/not/exist')

        self.assertRegex(str(e.exception),
                         r"^file '/does/not/exist' does not exist.")

        self.assertIn('QCJob.py', job._find_files(''))

        obs = job._system_call('ls tests/bin')
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
