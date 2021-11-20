import unittest
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import abspath, join
from functools import partial
from shutil import rmtree
import re


class TestJob(unittest.TestCase):
    def setUp(self):
        self.remove_these = []

    def tearDown(self):
        for some_path in self.remove_these:
            rmtree(some_path)

    def test_system_call(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

        output_dir = self.path('my_output_dir')

        job = Job(self.path('sample-sequence-directory'),
                  output_dir, '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], 1000, None)

        self.assertTrue(job._which('ls') in ['/bin/ls', '/usr/bin/ls'])

        with self.assertRaises(PipelineError) as e:
            job._file_check('/does/not/exist')

        self.assertRegex(str(e.exception),
                         r"^file '/does/not/exist' does not exist.")

        obs = job._find_files(package_root)
        obs = [re.sub(r'^.*?/sequence_processing_pipeline',
                      r'sequence_processing_pipeline', x) for x in obs]
        self.assertIn('sequence_processing_pipeline/QCJob.py', obs)

        obs = job._system_call('ls ' + join(package_root, 'tests', 'bin'))
        exp = {'stdout': 'bcl-convert\nbcl2fastq\nfastqc\n',
               'stderr': '',
               'return_code': 0}
        self.assertDictEqual(obs, exp)

        with self.assertRaises(PipelineError) as e:
            job._system_call('error')

        self.assertRegex(str(e.exception), (r"^Execute command-line statement"
                                            r" failure:\nCommand: error\nretur"
                                            r"n code: 127"))

        self.remove_these.append(output_dir)

    def test_group_commands(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

        job = Job(self.path('sample-sequence-directory'),
                  self.path('my_output_dir'), '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], 2, None)

        cmds = list(range(1, 8))
        cmds = [str(x) for x in cmds]
        results = job._group_commands(cmds)
        print(results)
        self.assertEqual(results[0], '1;3;5;7')
        self.assertEqual(results[1], '2;4;6')
        self.assertEqual(len(results), 2)


if __name__ == '__main__':
    unittest.main()
