import unittest
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import abspath, join, dirname
from os import makedirs
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

        job = Job(self.path('211021_A00000_0000_SAMPLE'),
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

        callback_results = []

        def my_callback(jid=None, status=None):
            callback_results.append((jid, status))

        obs = job._system_call('ls ' + join(package_root, 'tests', 'bin'),
                               callback=my_callback)

        exp = ['bcl2fastq\nbcl-convert\nfastqc\n',
               'bcl-convert\nbcl2fastq\nfastqc\n']

        self.assertIn(obs['stdout'], exp)
        self.assertEqual(obs['stderr'], '')
        self.assertEqual(obs['return_code'], 0)

        for item in callback_results:
            self.assertTrue(isinstance(item[0], int))
            self.assertIn(item[1], ['RUNNING', 'COMPLETED'])

        with self.assertRaises(PipelineError) as e:
            job._system_call('error')

        self.assertRegex(str(e.exception), (r'^Execute command-line statement'
                                            r' failure:\nCommand: '
                                            r'error\nreturn code:'
                                            r' 127'))

        self.remove_these.append(output_dir)

    def test_group_commands(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')

        job = Job(self.path('211021_A00000_0000_SAMPLE'),
                  self.path('my_output_dir'), '200nnn_xnnnnn_nnnn_xxxxxxxxxx',
                  ['ls'], 2, None)

        cmds = list(range(1, 8))
        cmds = [str(x) for x in cmds]
        results = job._group_commands(cmds)
        self.assertEqual(results[0], '1;3;5;7')
        self.assertEqual(results[1], '2;4;6')
        self.assertEqual(len(results), 2)

    def test_extract_project_names_from_fastq_dir(self):
        package_root = abspath('./sequence_processing_pipeline')
        base_path = partial(join, package_root, 'tests', 'data')

        dummy_fastqs = [
            "ConvertJob/NPH_15288/359180337_S27_L007_R1_001.fastq.gz",
            "ConvertJob/NPH_15288/BLANK2_NPH_2_E_S13_L007_R1_001.fastq.gz",
            "ConvertJob/Undetermined_S0_L007_R2_001.fastq.gz",
            ("NuQCJob/only-adapter-filtered/NPH_15288/"
             "359180398_S9_L007_R1_001.fastq.gz"),
            ("NuQCJob/only-adapter-filtered/NPH_15288/"
             "BLANK2_NPH_9_E_S69_L007_R2_001.fastq.gz"),
            ("NuQCJob/NPH_15288/filtered_sequences/"
             "359180396_S7_L007_R2_001.trimmed.fastq.gz")]

        self.remove_these.append(base_path('7b9d7d9c-2cd4-4d54-94ac-'
                                           '40e07a713585'))

        for fastq in dummy_fastqs:
            full_path = base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585',
                                  fastq)
            makedirs(dirname(full_path), exist_ok=True)
            with open(full_path, 'w') as f:
                f.write("This is a dummy file.")

        # fake a basic Job() object with enough metadata to test
        # extract_project_names_from_fastq_dir().
        job = Job(base_path('211021_A00000_0000_SAMPLE'),
                  base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585'),
                  '200nnn_xnnnnn_nnnn_xxxxxxxxxx', ['ls'], 2, None)

        # results from ConvertJob and NuQCJob should be equal.
        tmp = base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585', 'ConvertJob')
        obs = job.extract_project_names_from_fastq_dir(tmp)
        self.assertEqual(obs, ['NPH_15288'])

        tmp = base_path('7b9d7d9c-2cd4-4d54-94ac-40e07a713585', 'NuQCJob')
        obs = job.extract_project_names_from_fastq_dir(tmp)
        self.assertEqual(obs, ['NPH_15288'])


if __name__ == '__main__':
    unittest.main()
