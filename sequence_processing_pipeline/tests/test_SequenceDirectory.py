import unittest
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.PipelineError import PipelineError
import logging
from os.path import join
from os import getcwd


logging.basicConfig(level=logging.DEBUG)


class TestSequenceDirectory(unittest.TestCase):
    def test_job(self):
        pass