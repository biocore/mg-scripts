import unittest
from unittest.mock import patch
import tempfile
import gzip
import os
from os.path import join
from sequence_processing_pipeline.Commands import (split_similar_size_bins,
                                                   demux)
import io


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

        with tempfile.TemporaryDirectory() as tmp:
            exp = 2
            exp_1 = ('/foo/bar/a_R1_.fastq.gz\t/foo/bar/a_R2_.fastq.gz\tbar\n'
                     '/foo/bar/b_R1_.fastq.gz\t/foo/bar/b_R2_.fastq.gz\tbar\n')
            exp_2 = '/foo/baz/c_R1_.fastq.gz\t/foo/baz/c_R2_.fastq.gz\tbaz\n'

            stat.return_value = MockStat()  # 512MB
            glob.return_value = mockglob
            obs = split_similar_size_bins('foo', 1, tmp + '/prefix')
            self.assertEqual(obs, exp)
            obs_1 = open(tmp + '/prefix-1').read()
            self.assertEqual(obs_1, exp_1)
            obs_1 = open(tmp + '/prefix-2').read()
            self.assertEqual(obs_1, exp_2)

    def test_demux(self):
        with tempfile.TemporaryDirectory() as tmp:
            id_map = io.StringIO("1\ta_R1\ta_R2\tProject_12345\n"
                                 "2\tb_R1\tb_R2\tProject_12345\n")
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

            out_d = tmp
            encoded = '2'
            threads = 1

            demux(id_map, infile, out_d, encoded, threads)

            obs_r1 = gzip.open(join(tmp, 'Project_12345', 'b_R1.fastq.gz'),
                               'rt').read()
            obs_r2 = gzip.open(join(tmp, 'Project_12345', 'b_R2.fastq.gz'),
                               'rt').read()
            self.assertEqual(obs_r1, exp_data_r1)
            self.assertEqual(obs_r2, exp_data_r2)
            self.assertFalse(os.path.exists(join(tmp, 'a_R1.fastq.gz')))
            self.assertFalse(os.path.exists(join(tmp, 'a_R2.fastq.gz')))


if __name__ == '__main__':
    unittest.main()
