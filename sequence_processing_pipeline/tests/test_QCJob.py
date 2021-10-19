import shutil
import unittest
from os.path import exists, join, basename, abspath
from functools import partial
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.PipelineError import PipelineError
from os import makedirs
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
import logging

logging.basicConfig(level=logging.DEBUG)


class TestQCJob(unittest.TestCase):
    def setUp(self):
        # adjustable test_root helps w/testing in different environments.
        test_root = abspath('sequence_processing_pipeline')
        self.path = partial(join, test_root, 'tests', 'data')
        self.sample_sheet_path = self.path('good-sample-sheet.csv')
        self.mmi_db_path = self.path('mmi.db')
        self.sample_run_dir = self.path('MyRunDir')
        self.project_list = ['NYU_BMS_Melanoma_13059', 'Feist_11661',
                             'Gerwick_6123']
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.maxDiff = None

        try:
            shutil.rmtree(self.sample_run_dir)
        except FileNotFoundError:
            # sample_run_dir is created during testing. we don't mind if it
            # doesn't already exist, but any leftovers from previous runs
            # need to be removed.
            pass

        sheet = KLSampleSheet(self.sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        sample_ids = []
        for sample in valid_sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        self.sample_paths = []
        fastq_path = partial(join, self.sample_run_dir, 'Data', 'Fastq')
        for project_name in self.project_list:
            sample_path = fastq_path(project_name)
            makedirs(sample_path, exist_ok=True)
            self.sample_paths.append(sample_path)

            ids = [x[0] for x in
                   filter(lambda c: c[1] == project_name, sample_ids)]

            for id in ids:
                fr_fp = join(sample_path, f'{id}_R1_001.fastq.gz')
                rr_fp = join(sample_path, f'{id}_R2_001.fastq.gz')
                with open(fr_fp, 'w') as f:
                    f.write('This is a forward-read file.')
                with open(rr_fp, 'w') as f:
                    f.write('This is a reverse-read file.')

    def test_split_file_creation(self):
        for sample_path in self.sample_paths:
            qc_job = QCJob(self.sample_run_dir, self.sample_sheet_path,
                           self.mmi_db_path, 'queue_name', 1, 16, 24, '8gb',
                           'fastp', 'minimap2', 'samtools', [],
                           self.qiita_job_id, 30, sample_path)

        # assert that the Torque job files were created and are in the
        # proper location.
        exp = ['QCJob_1.sh', 'QCJob_2.sh', 'QCJob_3.sh']
        exp = set([join(self.sample_run_dir, x) for x in exp])
        obs = set([qc_job.script_paths[proj] for proj in qc_job.script_paths])
        self.assertEqual(obs, exp)

        # compare the expected content for QCJob_1-3.sh with the observed
        # results.
        remove_this = None
        for shell_file, lines_exp in [('QCJob_1.sh', self.exp_QCJob_1),
                                      ('QCJob_2.sh', self.exp_QCJob_2),
                                      ('QCJob_3.sh', self.exp_QCJob_3)]:
            obs_fp = [x for x in list(exp) if shell_file in x][0]
            with open(obs_fp, 'r') as f_obs:
                lines_obs = f_obs.readlines()
                lines_obs = [x.strip() for x in lines_obs]
                # Technically, remove_this only needs to be calculated
                # once. Assume there is only one line in the file that
                # begins w/cd.
                tmp = [x for x in lines_obs if x.startswith('cd ')][0]
                # remove the cd, but restore other spaces
                tmp = ' '.join(tmp.split(' ')[1:])
                # split off path prepending 'sequence_processing...'
                tmp = tmp.split('sequence_processing_pipeline/tests')
                # if tmp2 has '' at element zero, then nothing
                # prepends 'sequence_processing_pipeline' and nothing
                # needs to be removed.
                remove_this = tmp[0] if tmp[0] != '' else None
                if remove_this:
                    lines_obs = [x.replace(remove_this, '') for x in
                                 lines_obs]

                self.assertEqual(lines_obs, lines_exp)

        # assert that the array-details files were created and are in the
        # proper location.
        exp = ['split_file_Feist_11661.array-details',
               'split_file_Gerwick_6123.array-details']

        '''    'split_file_Feist_11661_0', 'split_file_Feist_11661_1',
               'split_file_Feist_11661_2', 'split_file_Feist_11661_3',
               'split_file_Gerwick_6123_0',
               'split_file_NYU_BMS_Melanoma_13059_0',
               'split_file_NYU_BMS_Melanoma_13059_1',
               'split_file_NYU_BMS_Melanoma_13059_2',
               'split_file_NYU_BMS_Melanoma_13059_3']
        '''

        exp = set([join(self.sample_run_dir, x) for x in exp])

        exp_map = {'split_file_Feist_11661.array-details':
                   (self.feist_11661_beg, self.feist_11661_end),
                   'split_file_Gerwick_6123.array-details':
                   (self.gerwic_6123_beg, self.gerwic_6123_end),
                   'split_file_NYU_BMS_Melanoma_13059.array-details':
                   (self.nyubms_mela_beg, self.nyubms_mela_end)}

        for some_path in exp:
            # assert files are in the proper location
            self.assertTrue(exists(some_path) is True)

            lines_exp_beg, lines_exp_end = exp_map[basename(some_path)]
            # compare the expected content for the paths in exp with the
            # observed results.
            with open(some_path, 'r') as f_obs:
                lines_obs = f_obs.readlines()
                lines_obs = [x.strip() for x in lines_obs]
                if remove_this:
                    lines_obs = [x.replace(remove_this, '') for x in lines_obs]
                lines_obs.sort()
                self.assertEqual(lines_obs[0], lines_exp_beg)
                self.assertEqual(lines_obs[-1], lines_exp_end)

    def test_qcjob_creation(self):
        with self.assertRaises(PipelineError) as e:
            # Pick the first project in the list and generate the path to the
            # Fastq directory.
            fastq_path = partial(join, self.sample_run_dir, 'Data', 'Fastq')
            sample_path = fastq_path(self.project_list[0])

            QCJob(self.sample_run_dir, 'not/path/to/sample/sheet',
                  self.mmi_db_path, 'queue_name', 1, 16, 24, '8gb', 'fastp',
                  'minimap2', 'samtools', [], self.qiita_job_id, 30,
                  sample_path)

        self.assertEqual(str(e.exception), "file 'not/path/to/sample/sheet' "
                                           "does not exist.")

    def test_generate_split_count(self):
        '''
        Extra testing for _generate_split_count() since the number of sample
        projects needed to exercise all of this method would be impractical.
        :return:
        '''
        fastq_path = partial(join, self.sample_run_dir, 'Data', 'Fastq')
        sample_path = fastq_path(self.project_list[0])

        qc_job = QCJob(self.sample_run_dir, self.sample_sheet_path,
                       self.mmi_db_path, 'queue_name', 1, 16, 24, '8gb',
                       'fastp', 'minimap2', 'samtools', [], self.qiita_job_id,
                       30, sample_path)

        input = [-1, 0, 1, 499, 500, 501, 999, 1000, 1001, 1999, 2000, 2001]
        exp_output = [1, 1, 1, 1, 1, 4, 4, 4, 10, 10, 10, 16]
        obs_output = [qc_job._generate_split_count(x) for x in input]

        self.assertEqual(obs_output, exp_output)

    exp_QCJob_1 = [
        '#!/bin/bash',
        ('#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_NYU_BMS_Melanoma_1305'
         '9'),
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        ('#PBS -o localhost:sequence_processing_pipeline/tests/data/MyRunDir/Q'
         'CJob_1.out.log.${PBS_ARRAYID}'),
        ('#PBS -e localhost:sequence_processing_pipeline/tests/data/MyRunDir/Q'
         'CJob_1.err.log.${PBS_ARRAYID}'),
        '#PBS -t 1-86%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/MyRunDir',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        ('cmd0=$(head -n $step sequence_processing_pipeline/tests/data/MyRunDi'
         'r/split_file_NYU_BMS_Melanoma_13059.array-details | tail -n 1)'),
        'eval $cmd0']

    exp_QCJob_2 = [
        '#!/bin/bash',
        '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_Feist_11661',
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        ('#PBS -o localhost:sequence_processing_pipeline/tests/data/MyRunDir/Q'
         'CJob_2.out.log.${PBS_ARRAYID}'),
        ('#PBS -e localhost:sequence_processing_pipeline/tests/data/MyRunDir/Q'
         'CJob_2.err.log.${PBS_ARRAYID}'),
        '#PBS -t 1-108%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/MyRunDir',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        ('cmd0=$(head -n $step sequence_processing_pipeline/tests/data/MyRunDi'
         'r/split_file_Feist_11661.array-details | tail -n 1)'),
        'eval $cmd0']

    exp_QCJob_3 = [
        '#!/bin/bash',
        '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_Gerwick_6123',
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        ('#PBS -o localhost:sequence_processing_pipeline/tests/data/MyRunDir/Q'
         'CJob_3.out.log.${PBS_ARRAYID}'),
        ('#PBS -e localhost:sequence_processing_pipeline/tests/data/MyRunDir/Q'
         'CJob_3.err.log.${PBS_ARRAYID}'),
        '#PBS -t 1-9%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/MyRunDir',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        ('cmd0=$(head -n $step sequence_processing_pipeline/tests/data/MyRunDi'
         'r/split_file_Gerwick_6123.array-details | tail -n 1)'),
        'eval $cmd0']

    feist_11661_beg = ('fastp --adapter_sequence AACC --adapter_sequence_r2 '
                       'GGTT -l 100 -i sequence_processing_pipeline/tests/da'
                       'ta/MyRunDir/Data/Fastq/Feist_11661/AB5075_AZM_TALE_i'
                       'n_MHB_A_baumannii_AB5075_WT_17_25_R1_001.fastq.gz -I'
                       ' sequence_processing_pipeline/tests/data/MyRunDir/Da'
                       'ta/Fastq/Feist_11661/AB5075_AZM_TALE_in_MHB_A_bauman'
                       'nii_AB5075_WT_17_25_R2_001.fastq.gz -w 16 -j sequenc'
                       'e_processing_pipeline/tests/data/MyRunDir/Data/Fastq'
                       '/Gerwick_6123/Feist_11661/Feist_11661/json/AB5075_AZ'
                       'M_TALE_in_MHB_A_baumannii_AB5075_WT_17_25_R1_001.jso'
                       'n -h sequence_processing_pipeline/tests/data/MyRunDi'
                       'r/Data/Fastq/Gerwick_6123/Feist_11661/Feist_11661/ht'
                       'ml/AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_2'
                       '5_R1_001.html -o sequence_processing_pipeline/tests/'
                       'data/MyRunDir/Data/Fastq/Gerwick_6123/Feist_11661/Fe'
                       'ist_11661/trimmed_sequences/AB5075_AZM_TALE_in_MHB_A'
                       '_baumannii_AB5075_WT_17_25_R1_001.fastp.fastq.gz -O '
                       'sequence_processing_pipeline/tests/data/MyRunDir/Dat'
                       'a/Fastq/Gerwick_6123/Feist_11661/Feist_11661/trimmed'
                       '_sequences/AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075'
                       '_WT_17_25_R2_001.fastp.fastq.gz -R AB5075_AZM_TALE_i'
                       'n_MHB_A_baumannii_AB5075_WT_17_25_R1_001_report')

    feist_11661_end = ('fastp --adapter_sequence AACC --adapter_sequence_r2 '
                       'GGTT -l 100 -i sequence_processing_pipeline/tests/da'
                       'ta/MyRunDir/Data/Fastq/Feist_11661/stALE_E_coli_A9_F'
                       '44_I1_R1_R1_001.fastq.gz -I sequence_processing_pipe'
                       'line/tests/data/MyRunDir/Data/Fastq/Feist_11661/stAL'
                       'E_E_coli_A9_F44_I1_R1_R2_001.fastq.gz -w 16 -j seque'
                       'nce_processing_pipeline/tests/data/MyRunDir/Data/Fas'
                       'tq/Gerwick_6123/Feist_11661/Feist_11661/json/stALE_E'
                       '_coli_A9_F44_I1_R1_R1_001.json -h sequence_processin'
                       'g_pipeline/tests/data/MyRunDir/Data/Fastq/Gerwick_61'
                       '23/Feist_11661/Feist_11661/html/stALE_E_coli_A9_F44_'
                       'I1_R1_R1_001.html -o sequence_processing_pipeline/te'
                       'sts/data/MyRunDir/Data/Fastq/Gerwick_6123/Feist_1166'
                       '1/Feist_11661/trimmed_sequences/stALE_E_coli_A9_F44_'
                       'I1_R1_R1_001.fastp.fastq.gz -O sequence_processing_p'
                       'ipeline/tests/data/MyRunDir/Data/Fastq/Gerwick_6123/'
                       'Feist_11661/Feist_11661/trimmed_sequences/stALE_E_co'
                       'li_A9_F44_I1_R1_R2_001.fastp.fastq.gz -R stALE_E_col'
                       'i_A9_F44_I1_R1_R1_001_report')

    gerwic_6123_beg = ('fastp --adapter_sequence AACC --adapter_sequence_r2 '
                       'GGTT -l 100 -i sequence_processing_pipeline/tests/da'
                       'ta/MyRunDir/Data/Fastq/Gerwick_6123/3A_R1_001.fastq.'
                       'gz -I sequence_processing_pipeline/tests/data/MyRunD'
                       'ir/Data/Fastq/Gerwick_6123/3A_R2_001.fastq.gz -w 16 '
                       '-j sequence_processing_pipeline/tests/data/MyRunDir/'
                       'Data/Fastq/Gerwick_6123/Gerwick_6123/Gerwick_6123/js'
                       'on/3A_R1_001.json -h sequence_processing_pipeline/te'
                       'sts/data/MyRunDir/Data/Fastq/Gerwick_6123/Gerwick_61'
                       '23/Gerwick_6123/html/3A_R1_001.html -o sequence_proc'
                       'essing_pipeline/tests/data/MyRunDir/Data/Fastq/Gerwi'
                       'ck_6123/Gerwick_6123/Gerwick_6123/trimmed_sequences/'
                       '3A_R1_001.fastp.fastq.gz -O sequence_processing_pipe'
                       'line/tests/data/MyRunDir/Data/Fastq/Gerwick_6123/Ger'
                       'wick_6123/Gerwick_6123/trimmed_sequences/3A_R2_001.f'
                       'astp.fastq.gz -R 3A_R1_001_report')

    gerwic_6123_end = ('fastp --adapter_sequence AACC --adapter_sequence_r2 '
                       'GGTT -l 100 -i sequence_processing_pipeline/tests/da'
                       'ta/MyRunDir/Data/Fastq/Gerwick_6123/ISB_R1_001.fastq'
                       '.gz -I sequence_processing_pipeline/tests/data/MyRun'
                       'Dir/Data/Fastq/Gerwick_6123/ISB_R2_001.fastq.gz -w 1'
                       '6 -j sequence_processing_pipeline/tests/data/MyRunDi'
                       'r/Data/Fastq/Gerwick_6123/Gerwick_6123/Gerwick_6123/'
                       'json/ISB_R1_001.json -h sequence_processing_pipeline'
                       '/tests/data/MyRunDir/Data/Fastq/Gerwick_6123/Gerwick'
                       '_6123/Gerwick_6123/html/ISB_R1_001.html -o sequence_'
                       'processing_pipeline/tests/data/MyRunDir/Data/Fastq/G'
                       'erwick_6123/Gerwick_6123/Gerwick_6123/trimmed_sequen'
                       'ces/ISB_R1_001.fastp.fastq.gz -O sequence_processing'
                       '_pipeline/tests/data/MyRunDir/Data/Fastq/Gerwick_612'
                       '3/Gerwick_6123/Gerwick_6123/trimmed_sequences/ISB_R2'
                       '_001.fastp.fastq.gz -R ISB_R1_001_report')

    nyubms_mela_beg = ('fastp --adapter_sequence AACC --adapter_sequence_r2 '
                       'GGTT -l 100 -i sequence_processing_pipeline/tests/da'
                       'ta/MyRunDir/Data/Fastq/NYU_BMS_Melanoma_13059/AP1733'
                       '01B04_R1_001.fastq.gz -I sequence_processing_pipelin'
                       'e/tests/data/MyRunDir/Data/Fastq/NYU_BMS_Melanoma_13'
                       '059/AP173301B04_R2_001.fastq.gz -w 16 -j sequence_pr'
                       'ocessing_pipeline/tests/data/MyRunDir/Data/Fastq/Ger'
                       'wick_6123/NYU_BMS_Melanoma_13059/NYU_BMS_Melanoma_13'
                       '059/json/AP173301B04_R1_001.json -h sequence_process'
                       'ing_pipeline/tests/data/MyRunDir/Data/Fastq/Gerwick_'
                       '6123/NYU_BMS_Melanoma_13059/NYU_BMS_Melanoma_13059/h'
                       'tml/AP173301B04_R1_001.html -o sequence_processing_p'
                       'ipeline/tests/data/MyRunDir/Data/Fastq/Gerwick_6123/'
                       'NYU_BMS_Melanoma_13059/NYU_BMS_Melanoma_13059/trimme'
                       'd_sequences/AP173301B04_R1_001.fastp.fastq.gz -O seq'
                       'uence_processing_pipeline/tests/data/MyRunDir/Data/F'
                       'astq/Gerwick_6123/NYU_BMS_Melanoma_13059/NYU_BMS_Mel'
                       'anoma_13059/trimmed_sequences/AP173301B04_R2_001.fas'
                       'tp.fastq.gz -R AP173301B04_R1_001_report')

    nyubms_mela_end = ('fastp --adapter_sequence AACC --adapter_sequence_r2 '
                       'GGTT -l 100 -i sequence_processing_pipeline/tests/da'
                       'ta/MyRunDir/Data/Fastq/NYU_BMS_Melanoma_13059/lp1278'
                       '96a01_R1_001.fastq.gz -I sequence_processing_pipelin'
                       'e/tests/data/MyRunDir/Data/Fastq/NYU_BMS_Melanoma_13'
                       '059/lp127896a01_R2_001.fastq.gz -w 16 -j sequence_pr'
                       'ocessing_pipeline/tests/data/MyRunDir/Data/Fastq/Ger'
                       'wick_6123/NYU_BMS_Melanoma_13059/NYU_BMS_Melanoma_13'
                       '059/json/lp127896a01_R1_001.json -h sequence_process'
                       'ing_pipeline/tests/data/MyRunDir/Data/Fastq/Gerwick_'
                       '6123/NYU_BMS_Melanoma_13059/NYU_BMS_Melanoma_13059/h'
                       'tml/lp127896a01_R1_001.html -o sequence_processing_p'
                       'ipeline/tests/data/MyRunDir/Data/Fastq/Gerwick_6123/'
                       'NYU_BMS_Melanoma_13059/NYU_BMS_Melanoma_13059/trimme'
                       'd_sequences/lp127896a01_R1_001.fastp.fastq.gz -O seq'
                       'uence_processing_pipeline/tests/data/MyRunDir/Data/F'
                       'astq/Gerwick_6123/NYU_BMS_Melanoma_13059/NYU_BMS_Mel'
                       'anoma_13059/trimmed_sequences/lp127896a01_R2_001.fas'
                       'tp.fastq.gz -R lp127896a01_R1_001_report')


if __name__ == '__main__':
    unittest.main()
