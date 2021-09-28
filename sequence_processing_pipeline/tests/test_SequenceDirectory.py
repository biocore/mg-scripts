import unittest
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.PipelineError import PipelineError
import os
import logging


logging.basicConfig(level=logging.DEBUG)


class TestSequenceDirectory(unittest.TestCase):
    def test_job(self):
        # Job is a minimal base-class. Use lsJob to test basic Job
        # functionality.

        # base_path = join('sequence_processing_pipeline', 'tests', 'data')

        # raise PipelineError if directory path does not exist.
        self.assertRaises(PipelineError, SequenceDirectory,
                          'tests/data/some_other_directory',
                          'tests/data/good-sample-sheet.csv')
        # raise PipelineError if sample file does not exist.
        self.assertRaises(PipelineError, SequenceDirectory,
                          'tests/data/sample-sequence-directory',
                          'tests/data/file_does_not_exist.txt')
        # raise PipelineError if both don't exist.
        self.assertRaises(PipelineError, SequenceDirectory,
                          'tests/data/some_other_directory',
                          'tests/data/file_does_not_exist.txt')

        # verify 'Data/Fastq' subdirectory is created.
        sdo = SequenceDirectory('tests/data/sample-sequence-directory',
                                'tests/data/good-sample-sheet.csv')
        s = os.path.exists('tests/data/sample-sequence-directory/Data/Fastq')
        self.assertTrue(s)

        # assert should not raise for directory 'Data/Fastq' already existing
        # when creating a new object.
        success = True
        try:
            another = SequenceDirectory('tests/data/sample-sequence-directory',
                                        'tests/data/good-sample-sheet.csv')
        except PipelineError:
            success = False

        self.assertTrue(success)

        # assert that SequenceDirectory object correctly extracts the metadata
        # from a known sample sheet.

        results = sdo.process()
        self.assertEqual(results['chemistry'], 'Default')
        self.assertEqual(results['base_mask'],
                         '--use-bases-mask Y150,I12,I12,Y150')
        self.assertEqual(results['experiment_name'], 'RKL0042')
        self.assertEqual(results['sequence_directory'],
                         'tests/data/sample-sequence-directory')
        self.assertEqual(results['sample_sheet_path'],
                         'tests/data/good-sample-sheet.csv')

        # assert that SequenceDirectory object correctly raises an error when
        # initialized with an invalid sample sheet.
        a = 'tests/data/sample-sequence-directory'
        b = 'tests/data/no-project-name-sample-sheet.csv'
        another = SequenceDirectory(a, b)

        self.assertRaises(PipelineError, another.process)
