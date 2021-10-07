from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import logging
import unittest
from os import makedirs
from os.path import join, dirname, abspath
from functools import partial


logging.basicConfig(level=logging.DEBUG)


class TestPipeline(unittest.TestCase):
    def test_creation(self):
        pass