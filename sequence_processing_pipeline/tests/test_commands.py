import unittest
from unittest.mock import patch
from tempfile import TemporaryDirectory
import gzip
import os
from sequence_processing_pipeline.Commands import (split_similar_size_bins,
                                                   demux)
import io
from os.path import join


class CommandTests(unittest.TestCase):
    @patch('os.stat')
    @patch('glob.glob')
    def test_split_similar_size_bins(self, glob, stat):
        class MockStat:
            st_size = 2 ** 28  # 256MB

        mockglob = ['/foo/bar/a_R1_.fastq.gz',
                    '/foo/bar/b_R2_.fastq.gz',
                    '/foo/bar/a_R2_.fastq.gz',
                    '/foo/baz/c_R2_.fastq.gz',
                    '/foo/baz/c_R1_.fastq.gz',
                    '/foo/bar/b_R1_.fastq.gz']

        with TemporaryDirectory() as tmp:
            exp = (2, 1073741824)
            stat.return_value = MockStat()  # 512MB
            glob.return_value = mockglob
            obs = split_similar_size_bins('foo', 1, tmp + '/prefix')
            self.assertEqual(obs, exp)

            exp_1 = ('/foo/bar/a_R1_.fastq.gz\t/foo/bar/a_R2_.fastq.gz\tbar\n'
                     '/foo/bar/b_R1_.fastq.gz\t/foo/bar/b_R2_.fastq.gz\tbar\n')
            exp_2 = '/foo/baz/c_R1_.fastq.gz\t/foo/baz/c_R2_.fastq.gz\tbaz\n'

            obs_1 = open(tmp + '/prefix-1').read()
            self.assertEqual(obs_1, exp_1)
            obs_1 = open(tmp + '/prefix-2').read()
            self.assertEqual(obs_1, exp_2)

    def test_demux(self):
        with TemporaryDirectory() as tmp:
            id_map = [
                        ["1", "a_R1", "a_R2", "Project_12345"],
                        ["2", "b_R1", "b_R2", "Project_12345"]
                        ]

            infile_data = '\n'.join(['@1::MUX::foo/1', 'ATGC', '+', '!!!!',
                                     '@1::MUX::foo/2', 'ATGC', '+', '!!!!',
                                     '@1::MUX::bar/1', 'ATGC', '+', '!!!!',
                                     '@1::MUX::bar/2', 'ATGC', '+', '!!!!',
                                     '@2::MUX::baz/1', 'ATGC', '+', '!!!!',
                                     '@2::MUX::baz/2', 'ATGC', '+', '!!!!',
                                     '@2::MUX::bing/1', 'ATGC', '+', '!!!!',
                                     '@2::MUX::bing/2', 'ATGC', '+', '!!!!',
                                     ''])
            infile = io.StringIO(infile_data)
            exp_data_r1 = '\n'.join(['@baz/1', 'ATGC', '+', '!!!!',
                                     '@bing/1', 'ATGC', '+', '!!!!', ''])
            exp_data_r2 = '\n'.join(['@baz/2', 'ATGC', '+', '!!!!',
                                     '@bing/2', 'ATGC', '+', '!!!!', ''])

            exp_data_r1 = ['@baz/1', 'ATGC', '+', '!!!!',
                           '@bing/1', 'ATGC', '+', '!!!!']
            exp_data_r2 = ['@baz/2', 'ATGC', '+', '!!!!',
                           '@bing/2', 'ATGC', '+', '!!!!']

            task = 1
            maxtask = 2

            demux(id_map, infile, tmp, task, maxtask)

            obs_r1 = gzip.open(join(tmp, 'Project_12345', 'b_R1.fastq.gz'),
                               'rt').read()
            obs_r2 = gzip.open(join(tmp, 'Project_12345', 'b_R2.fastq.gz'),
                               'rt').read()
            exp = '\n'.join(exp_data_r1) + '\n'
            self.assertEqual(obs_r1, exp)

            exp = '\n'.join(exp_data_r2) + '\n'
            self.assertEqual(obs_r2, exp)

            self.assertFalse(os.path.exists(join(tmp, 'a_R1.fastq.gz')))
            self.assertFalse(os.path.exists(join(tmp, 'a_R2.fastq.gz')))


if __name__ == '__main__':
    unittest.main()
