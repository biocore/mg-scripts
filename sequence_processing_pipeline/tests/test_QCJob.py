import shutil
import unittest
from os.path import join, abspath, exists, basename
from functools import partial
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.PipelineError import PipelineError
from os import makedirs
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet


class TestQCJob(unittest.TestCase):
    def setUp(self):
        # adjustable test_root helps w/testing in different environments.
        package_root = abspath('sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')
        self.sample_sheet_path = self.path('good-sample-sheet.csv')
        self.mmi_db_path = self.path('mmi.db')
        self.project_list = ['NYU_BMS_Melanoma_13059', 'Feist_11661',
                             'Gerwick_6123']
        self.qiita_job_id = 'abcdabcdabcdabcdabcdabcdabcdabcd'
        self.maxDiff = None
        self.output_path = self.path('output_dir')
        self.fastq_root_path = join(self.output_path, 'ConvertJob')

        try:
            shutil.rmtree(self.output_path)
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

        self.fastq_path = partial(join, self.output_path, 'ConvertJob')
        for project_name in self.project_list:
            # strip the qiita-id from a project-name in order to test
            # QCJob's ability to match project directories for both new and
            # legacy project folders w/in a run-directory.
            tmp = 'Feist' if project_name == 'Feist_11661' else project_name
            sample_path = self.fastq_path(tmp)
            makedirs(sample_path, exist_ok=True)

            ids = [x[0] for x in
                   filter(lambda c: c[1] == project_name, sample_ids)]

            for id in ids:
                fr_fp = join(sample_path, f'{id}_R1_001.fastq.gz')
                rr_fp = join(sample_path, f'{id}_R2_001.fastq.gz')
                with open(fr_fp, 'w') as f:
                    f.write('This is a forward-read file.')
                with open(rr_fp, 'w') as f:
                    f.write('This is a reverse-read file.')

    def tearDown(self):
        shutil.rmtree(self.output_path)

    def test_qcjob_creation(self):
        with self.assertRaises(PipelineError) as e:
            QCJob(self.fastq_root_path, self.output_path,
                  'not/path/to/sample/sheet', self.mmi_db_path,
                  'queue_name', 1, 16, 24, '8gb', 'fastp', 'minimap2',
                  'samtools', [], self.qiita_job_id, 30, 1000)

        self.assertEqual(str(e.exception), "file 'not/path/to/sample/sheet' "
                                           "does not exist.")

    def test_split_file_creation(self):
        qc_job = QCJob(self.fastq_root_path, self.output_path,
                       self.sample_sheet_path, self.mmi_db_path, 'queue_name',
                       1, 16, 24, '8gb', 'fastp', 'minimap2', 'samtools', [],
                       self.qiita_job_id, 30, 1000)

        # assert that the Torque job files were created and are in the
        # proper location.
        exp = ['QCJob_Feist_11661.sh', 'QCJob_Gerwick_6123.sh',
               'QCJob_NYU_BMS_Melanoma_13059.sh']
        exp = set([join(qc_job.output_path, x) for x in exp])
        obs = set([qc_job.script_paths[proj] for proj in qc_job.script_paths])
        self.assertEqual(obs, exp)

        # compare the expected content for QCJob_1-3.sh with the observed
        # results.
        remove_this = None
        for shell_file, lines_exp in [('QCJob_NYU_BMS_Melanoma_13059.sh',
                                       self.exp_QCJob_1),
                                      ('QCJob_Feist_11661.sh',
                                       self.exp_QCJob_2),
                                      ('QCJob_Gerwick_6123.sh',
                                       self.exp_QCJob_3)]:
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
        exp = ['QCJob_Feist_11661.array-details',
               'QCJob_Gerwick_6123.array-details',
               'QCJob_NYU_BMS_Melanoma_13059.array-details']

        exp = set([join(qc_job.output_path, x) for x in exp])

        for some_path in exp:
            # assert files are in the proper location
            self.assertTrue(exists(some_path) is True)

            config = self.exp_map[basename(some_path)]
            # compare the expected content for the paths in exp with the
            # observed results.
            with open(some_path, 'r') as f_obs:
                lines_obs = f_obs.readlines()
                lines_obs = [x.strip() for x in lines_obs]
                if remove_this:
                    lines_obs = [x.replace(remove_this, '') for x in lines_obs]
                lines_obs.sort()
                self.assertEqual(lines_obs[0], config['first_line'])
                self.assertEqual(lines_obs[-1], config['last_line'])
                self.assertEqual(len(lines_obs), config['count'])

    def test_command_chaining(self):
        # limit job array length to 250 jobs to test chaining of multiple
        # commands as a single command.

        qc_job = QCJob(self.fastq_root_path, self.output_path,
                       self.sample_sheet_path, self.mmi_db_path, 'queue_name',
                       1, 16, 24, '8gb', 'fastp', 'minimap2', 'samtools', [],
                       self.qiita_job_id, 30, 250)

        # assert that the Torque job files were created and are in the
        # proper location.
        exp = ['QCJob_Feist_11661.sh', 'QCJob_Gerwick_6123.sh',
               'QCJob_NYU_BMS_Melanoma_13059.sh']
        exp = set([join(qc_job.output_path, x) for x in exp])
        obs = set([qc_job.script_paths[proj] for proj in qc_job.script_paths])
        self.assertEqual(obs, exp)

        # compare the expected content for QCJob_1-3.sh with the observed
        # results.
        remove_this = None
        for shell_file, lines_exp in [('QCJob_NYU_BMS_Melanoma_13059.sh',
                                       self.exp_QCJob_1_chained),
                                      ('QCJob_Feist_11661.sh',
                                       self.exp_QCJob_2_chained),
                                      ('QCJob_Gerwick_6123.sh',
                                       self.exp_QCJob_3_chained)]:
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
        exp = ['QCJob_Feist_11661.array-details',
               'QCJob_Gerwick_6123.array-details',
               'QCJob_NYU_BMS_Melanoma_13059.array-details']

        exp = set([join(qc_job.output_path, x) for x in exp])

        for some_path in exp:
            # assert files are in the proper location
            self.assertTrue(exists(some_path) is True)

            config = self.exp_map_chained[basename(some_path)]
            # compare the expected content for the paths in exp with the
            # observed results.
            with open(some_path, 'r') as f_obs:
                lines_obs = f_obs.readlines()
                lines_obs = [x.strip() for x in lines_obs]
                if remove_this:
                    lines_obs = [x.replace(remove_this, '') for x in lines_obs]
                    if lines_obs[0] != config['first_line']:
                        print("FIRST LINE OBSERVED ##################")
                        print(lines_obs[0])
                        print("FIRST LINE EXPECTED ##################")
                        print(config['first_line'])
                        print("##################")
                    if lines_obs[-1] != config['last_line']:
                        print("LAST LINE OBSERVED ##################")
                        print(lines_obs[-1])
                        print("LAST LINE EXPECTED ##################")
                        print(config['last_line'])
                        print("##################")
                self.assertEqual(lines_obs[0], config['first_line'])
                self.assertEqual(lines_obs[-1], config['last_line'])
                self.assertEqual(len(lines_obs), config['count'])

    exp_QCJob_1 = [
        '#!/bin/bash',
        ('#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_NYU_BMS_Melanoma_1305'
         '9'),
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        '#PBS -t 1-384%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/output_dir/QCJob',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        'cmd0=$(head -n $step sequence_processing_pipeline/tests/data/output_d'
        'ir/QCJob/QCJob_NYU_BMS_Melanoma_13059.array-details | tail -n 1)',
        'eval $cmd0']

    exp_QCJob_2 = [
        '#!/bin/bash',
        '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_Feist_11661',
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        '#PBS -t 1-432%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/output_dir/QCJob',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        'cmd0=$(head -n $step sequence_processing_pipeline/tests/data/output_d'
        'ir/QCJob/QCJob_Feist_11661.array-details | tail -n 1)',
        'eval $cmd0']

    exp_QCJob_3 = [
        '#!/bin/bash',
        '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_Gerwick_6123',
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        '#PBS -t 1-9%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/output_dir/QCJob',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        'cmd0=$(head -n $step sequence_processing_pipeline/tests/data/output_d'
        'ir/QCJob/QCJob_Gerwick_6123.array-details | tail -n 1)',
        'eval $cmd0']

    exp_map = {'QCJob_Feist_11661.array-details': {
                   'first_line': "fastp --adapter_sequence AACC "
                                 "--adapter_sequence_r2 GGTT -l 100 -i sequenc"
                                 "e_processing_pipeline/tests/data/output_dir/"
                                 "ConvertJob/Feist/AB5075_AZM_TALE_in_MH"
                                 "B_A_baumannii_AB5075_WT_17_25_R1_001.fastq.g"
                                 "z -I sequence_processing_pipeline/tests/data"
                                 "/output_dir/ConvertJob/Feist/AB5075_AZ"
                                 "M_TALE_in_MHB_A_baumannii_AB5075_WT_17_25_R2"
                                 "_001.fastq.gz -w 16 -j sequence_processing_p"
                                 "ipeline/tests/data/output_dir/QCJob"
                                 "/Feist_11661/json/AB5075_AZM_TALE_in_MHB_"
                                 "A_baumannii_AB5075_WT_17_25_R1_001.json -h s"
                                 "equence_processing_pipeline/tests/data/outpu"
                                 "t_dir/QCJob/Feist_11661/html/AB5"
                                 "075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17"
                                 "_25_R1_001.html -o sequence_processing_pipel"
                                 "ine/tests/data/output_dir/QCJob/"
                                 "Feist_11661/trimmed_sequences/AB5075_AZM_TAL"
                                 "E_in_MHB_A_baumannii_AB5075_WT_17_25_R1_001."
                                 "fastp.fastq.gz -O sequence_processing_pipeli"
                                 "ne/tests/data/output_dir/QCJob/F"
                                 "eist_11661/trimmed_sequences/AB5075_AZM_TALE"
                                 "_in_MHB_A_baumannii_AB5075_WT_17_25_R2_001.f"
                                 "astp.fastq.gz -R AB5075_AZM_TALE_in_MHB_A_ba"
                                 "umannii_AB5075_WT_17_25_R1_001_report",
                   'last_line':  "fastp --adapter_sequence AACC "
                                 "--adapter_sequence_r2 GGTT -l 100 -i sequenc"
                                 "e_processing_pipeline/tests/data/output_dir/"
                                 "ConvertJob/Feist/stALE_E_coli_A9_F44_I"
                                 "1_R1_R2_001.fastq.gz -I sequence_processing_"
                                 "pipeline/tests/data/output_dir/ConvertJob/Fe"
                                 "ist/stALE_E_coli_A9_F44_I1_R1_R2_001.f"
                                 "astq.gz -w 16 -j sequence_processing_pipelin"
                                 "e/tests/data/output_dir/QCJob/Fe"
                                 "ist_11661/json/stALE_E_coli_A9_F44_I1_R1_R2_"
                                 "001.json -h sequence_processing_pipeline/tes"
                                 "ts/data/output_dir/QCJob/Feist_1"
                                 "1661/html/stALE_E_coli_A9_F44_I1_R1_R2_001.h"
                                 "tml -o sequence_processing_pipeline/tests/da"
                                 "ta/output_dir/QCJob/Feist_11661/"
                                 "trimmed_sequences/stALE_E_coli_A9_F44_I1_R1_"
                                 "R2_001.fastp.fastq.gz -O sequence_processing"
                                 "_pipeline/tests/data/output_dir/QCJob/Feist_"
                                 "11661/trimmed_sequences/stALE_E_"
                                 "coli_A9_F44_I1_R1_R2_001.fastp.fastq.gz -R "
                                 "stALE_E_coli_A9_F44_I1_R1_R2_001_report",
                   'count': 432},

               'QCJob_Gerwick_6123.array-details': {
                   'first_line': "fastp --adapter_sequence AACC --adapter_sequ"
                                 "ence_r2 GGTT -l 100 -i sequence_processing_p"
                                 "ipeline/tests/data/output_dir/ConvertJob/Ger"
                                 "wick_6123/3A_R1_001.fastq.gz -I sequence_pro"
                                 "cessing_pipeline/tests/data/output_dir/Conve"
                                 "rtJob/Gerwick_6123/3A_R2_001.fastq.gz -w 16 "
                                 "-j sequence_processing_pipeline/tests/data/o"
                                 "utput_dir/QCJob/Gerwick_6123/fastp_reports_d"
                                 "ir/json/3A_R1_001.json -h sequence_processin"
                                 "g_pipeline/tests/data/output_dir/QCJob/Gerwi"
                                 "ck_6123/fastp_reports_dir/html/3A_R1_001.htm"
                                 "l --stdout | minimap2 -ax sr -t 16 sequence_"
                                 "processing_pipeline/tests/data/mmi.db - -a |"
                                 " samtools fastq -@ 16 -f 12 -F 256 -1 sequen"
                                 "ce_processing_pipeline/tests/data/output_dir"
                                 "/QCJob/Gerwick_6123/filtered_sequences/3A_R1"
                                 "_001.trimmed.fastq.gz -2 sequence_processing"
                                 "_pipeline/tests/data/output_dir/QCJob/Gerwic"
                                 "k_6123/filtered_sequences/3A_R2_001.trimmed."
                                 "fastq.gz",

                   'last_line': "fastp --adapter_sequence AACC --adapter_seque"
                                "nce_r2 GGTT -l 100 -i sequence_processing_pip"
                                "eline/tests/data/output_dir/ConvertJob/Gerwic"
                                "k_6123/ISB_R1_001.fastq.gz -I sequence_proces"
                                "sing_pipeline/tests/data/output_dir/ConvertJo"
                                "b/Gerwick_6123/ISB_R2_001.fastq.gz -w 16 -j s"
                                "equence_processing_pipeline/tests/data/output"
                                "_dir/QCJob/Gerwick_6123/fastp_reports_dir/jso"
                                "n/ISB_R1_001.json -h sequence_processing_pipe"
                                "line/tests/data/output_dir/QCJob/Gerwick_6123"
                                "/fastp_reports_dir/html/ISB_R1_001.html --std"
                                "out | minimap2 -ax sr -t 16 sequence_processi"
                                "ng_pipeline/tests/data/mmi.db - -a | samtools"
                                " fastq -@ 16 -f 12 -F 256 -1 sequence_process"
                                "ing_pipeline/tests/data/output_dir/QCJob/Gerw"
                                "ick_6123/filtered_sequences/ISB_R1_001.trimme"
                                "d.fastq.gz -2 sequence_processing_pipeline/te"
                                "sts/data/output_dir/QCJob/Gerwick_6123/filter"
                                "ed_sequences/ISB_R2_001.trimmed.fastq.gz",
                   'count': 9},
               'QCJob_NYU_BMS_Melanoma_13059.array-details': {
                   'first_line': "fastp --adapter_sequence AACC "
                                 "--adapter_sequence_r2 GGTT -l 100 -i sequenc"
                                 "e_processing_pipeline/tests/data/output_dir/"
                                 "ConvertJob/NYU_BMS_Melanoma_13059/22_001_710"
                                 "_503_791_00_R1_001.fastq.gz -I sequence_proc"
                                 "essing_pipeline/tests/data/output_dir/Conver"
                                 "tJob/NYU_BMS_Melanoma_13059/22_001_710_503_7"
                                 "91_00_R2_001.fastq.gz -w 16 -j sequence_proc"
                                 "essing_pipeline/tests/data/output_dir/QCJob/"
                                 "NYU_BMS_Melanoma_1305"
                                 "9/json/22_001_710_503_791_00_R1_001.json -h "
                                 "sequence_processing_pipeline/tests/data/outp"
                                 "ut_dir/QCJob/NYU_BMS_"
                                 "Melanoma_13059/html/22_001_710_503_791_00_R1"
                                 "_001.html -o sequence_processing_pipeline/te"
                                 "sts/data/output_dir/QCJob/NYU_BMS_Melanoma_1"
                                 "3059/trimmed_sequence"
                                 "s/22_001_710_503_791_00_R1_001.fastp.fastq.g"
                                 "z -O sequence_processing_pipeline/tests/data"
                                 "/output_dir/QCJob/NYU"
                                 "_BMS_Melanoma_13059/trimmed_sequences/22_001"
                                 "_710_503_791_00_R2_001.fastp.fastq.gz -R "
                                 "22_001_710_503_791_00_R1_001_report",
                   'last_line': "fastp --adapter_sequence AACC "
                                 "--adapter_sequence_r2 GGTT -l 100 -i sequenc"
                                 "e_processing_pipeline/tests/data/output_dir/"
                                 "ConvertJob/NYU_BMS_Melanoma_13059/lp127896a0"
                                 "1_R1_001.fastq.gz -I sequence_processing_pip"
                                 "eline/tests/data/output_dir/ConvertJob/NYU_B"
                                 "MS_Melanoma_13059/lp127896a01_R2_001.fastq.g"
                                 "z -w 16 -j sequence_processing_pipeline/test"
                                 "s/data/output_dir/QCJob/NYU_BMS_Melanoma_130"
                                 "59/json/lp127896a01_R"
                                 "1_001.json -h sequence_processing_pipeline/t"
                                 "ests/data/output_dir/QCJob/NYU_BMS_Melanoma_"
                                 "13059/html/lp127896a0"
                                 "1_R1_001.html -o sequence_processing_pipelin"
                                 "e/tests/data/output_dir/QCJob/NYU_BMS_Melano"
                                 "ma_13059/trimmed_sequ"
                                 "ences/lp127896a01_R1_001.fastp.fastq.gz -O s"
                                 "equence_processing_pipeline/tests/data/outpu"
                                 "t_dir/QCJob/NYU_BMS_M"
                                 "elanoma_13059/trimmed_sequences/lp127896a01_"
                                 "R2_001.fastp.fastq.gz -R "
                                 "lp127896a01_R1_001_report",
                   'count': 384}
               }

    exp_QCJob_1_chained = [
        '#!/bin/bash',
        ('#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_NYU_BMS_Melanoma_1305'
         '9'),
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        '#PBS -t 1-250%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/output_dir/QCJob',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        'cmd0=$(head -n $step sequence_processing_pipeline/tests/data/output_d'
        'ir/QCJob/QCJob_NYU_BMS_Melanoma_13059.array-details | tail -n 1)',
        'eval $cmd0']

    exp_QCJob_2_chained = [
        '#!/bin/bash',
        '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_Feist_11661',
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        '#PBS -t 1-250%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/output_dir/QCJob',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        'cmd0=$(head -n $step sequence_processing_pipeline/tests/data/output_d'
        'ir/QCJob/QCJob_Feist_11661.array-details | tail -n 1)',
        'eval $cmd0']

    exp_QCJob_3_chained = [
        '#!/bin/bash',
        '#PBS -N abcdabcdabcdabcdabcdabcdabcdabcd_QCJob_Gerwick_6123',
        '#PBS -q queue_name',
        '#PBS -l nodes=1:ppn=16',
        '#PBS -V',
        '#PBS -l walltime=24:00:00',
        '#PBS -l mem=8gb',
        '#PBS -t 1-9%30',
        'set -x',
        'date',
        'hostname',
        'echo ${PBS_JOBID} ${PBS_ARRAYID}',
        'cd sequence_processing_pipeline/tests/data/output_dir/QCJob',
        'offset=${PBS_ARRAYID}',
        'step=$(( $offset - 0 ))',
        'cmd0=$(head -n $step sequence_processing_pipeline/tests/data/output_d'
        'ir/QCJob/QCJob_Gerwick_6123.array-details | tail -n 1)',
        'eval $cmd0']

    exp_map_chained = {'QCJob_Feist_11661.array-details': {
        'first_line': "fastp --adapter_sequence AACC --adapter_sequence_r2 GGT"
                      "T -l 100 -i sequence_processing_pipeline/tests/data/out"
                      "put_dir/ConvertJob/Feist/AB5075_AZM_TALE_in_MHB_A_bauma"
                      "nnii_AB5075_WT_17_25_R1_001.fastq.gz -I sequence_proces"
                      "sing_pipeline/tests/data/output_dir/ConvertJob/Feist/AB"
                      "5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25_R2_001"
                      ".fastq.gz -w 16 -j sequence_processing_pipeline/tests/d"
                      "ata/output_dir/QCJob/Feist_11661/json/AB507"
                      "5_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25_R1_001.js"
                      "on -h sequence_processing_pipeline/tests/data/output_di"
                      "r/QCJob/Feist_11661/html/AB5075_AZM_TALE_in"
                      "_MHB_A_baumannii_AB5075_WT_17_25_R1_001.html -o sequenc"
                      "e_processing_pipeline/tests/data/output_dir/QCJob/Feist"
                      "_11661/trimmed_sequences/AB5075_AZM_TALE_in"
                      "_MHB_A_baumannii_AB5075_WT_17_25_R1_001.fastp.fastq.gz "
                      "-O sequence_processing_pipeline/tests/data/output_dir/Q"
                      "CJob/Feist_11661/trimmed_sequences/AB5075_A"
                      "ZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25_R2_001.fastp"
                      ".fastq.gz -R AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_"
                      "WT_17_25_R1_001_report;fastp --adapter_sequence AACC --"
                      "adapter_sequence_r2 GGTT -l 100 -i sequence_processing_"
                      "pipeline/tests/data/output_dir/ConvertJob/Feist/P21_E_c"
                      "oli_ELI344_R1_001.fastq.gz -I sequence_processing_pipel"
                      "ine/tests/data/output_dir/ConvertJob/Feist/P21_E_coli_E"
                      "LI344_R2_001.fastq.gz -w 16 -j sequence_processing_pipe"
                      "line/tests/data/output_dir/QCJob/Feist_1166"
                      "1/json/P21_E_coli_ELI344_R1_001.json -h sequence_proces"
                      "sing_pipeline/tests/data/output_dir/QCJob/F"
                      "eist_11661/html/P21_E_coli_ELI344_R1_001.html -o sequen"
                      "ce_processing_pipeline/tests/data/output_dir/QCJob/Feis"
                      "t_11661/trimmed_sequences/P21_E_coli_ELI344"
                      "_R1_001.fastp.fastq.gz -O sequence_processing_pipeline/"
                      "tests/data/output_dir/QCJob/Feist_11661/tri"
                      "mmed_sequences/P21_E_coli_ELI344_R2_001.fastp.fastq.gz "
                      "-R P21_E_coli_ELI344_R1_001_report",
        'last_line': "fastp --adapter_sequence AACC --adapter_sequence_r2 GGTT"
                     " -l 100 -i sequence_processing_pipeline/tests/data/outpu"
                     "t_dir/ConvertJob/Feist/JM-Metabolic__GN0_2404_R1_001.fas"
                     "tq.gz -I sequence_processing_pipeline/tests/data/output_"
                     "dir/ConvertJob/Feist/JM-Metabolic__GN0_2404_R2_001.fastq"
                     ".gz -w 16 -j sequence_processing_pipeline/tests/data/out"
                     "put_dir/QCJob/Feist_11661/json/JM-Metabolic_"
                     "_GN0_2404_R1_001.json -h sequence_processing_pipeline/te"
                     "sts/data/output_dir/QCJob/Feist_11661/html/J"
                     "M-Metabolic__GN0_2404_R1_001.html -o sequence_processing"
                     "_pipeline/tests/data/output_dir/QCJob/Feist_"
                     "11661/trimmed_sequences/JM-Metabolic__GN0_2404_R1_001.fa"
                     "stp.fastq.gz -O sequence_processing_pipeline/tests/data/"
                     "output_dir/QCJob/Feist_11661/trimmed_sequenc"
                     "es/JM-Metabolic__GN0_2404_R2_001.fastp.fastq.gz -R JM-Me"
                     "tabolic__GN0_2404_R1_001_report",
        'count': 250},
        'QCJob_Gerwick_6123.array-details': {
            'first_line': "fastp --adapter_sequence AACC --adapter_sequence_r2"
                          " GGTT -l 100 -i sequence_processing_pipeline/tests/"
                          "data/output_dir/ConvertJob/Gerwick_6123/3A_R1_001.f"
                          "astq.gz -I sequence_processing_pipeline/tests/data/"
                          "output_dir/ConvertJob/Gerwick_6123/3A_R2_001.fastq."
                          "gz -w 16 -j sequence_processing_pipeline/tests/data"
                          "/output_dir/QCJob/Gerwick_6123/fastp_reports_dir/js"
                          "on/3A_R1_001.json -h sequence_processing_pipeline/t"
                          "ests/data/output_dir/QCJob/Gerwick_6123/fastp_repor"
                          "ts_dir/html/3A_R1_001.html --stdout | minimap2 -ax "
                          "sr -t 16 sequence_processing_pipeline/tests/data/mm"
                          "i.db - -a | samtools fastq -@ 16 -f 12 -F 256 -1 se"
                          "quence_processing_pipeline/tests/data/output_dir/QC"
                          "Job/Gerwick_6123/filtered_sequences/3A_R1_001.trimm"
                          "ed.fastq.gz -2 sequence_processing_pipeline/tests/d"
                          "ata/output_dir/QCJob/Gerwick_6123/filtered_sequence"
                          "s/3A_R2_001.trimmed.fastq.gz",
            'last_line': "fastp --adapter_sequence AACC --adapter_sequence_r2 "
                         "GGTT -l 100 -i sequence_processing_pipeline/tests/da"
                         "ta/output_dir/ConvertJob/Gerwick_6123/ISB_R1_001.fas"
                         "tq.gz -I sequence_processing_pipeline/tests/data/out"
                         "put_dir/ConvertJob/Gerwick_6123/ISB_R2_001.fastq.gz "
                         "-w 16 -j sequence_processing_pipeline/tests/data/out"
                         "put_dir/QCJob/Gerwick_6123/fastp_reports_dir/json/IS"
                         "B_R1_001.json -h sequence_processing_pipeline/tests/"
                         "data/output_dir/QCJob/Gerwick_6123/fastp_reports_dir"
                         "/html/ISB_R1_001.html --stdout | minimap2 -ax sr -t "
                         "16 sequence_processing_pipeline/tests/data/mmi.db - "
                         "-a | samtools fastq -@ 16 -f 12 -F 256 -1 sequence_p"
                         "rocessing_pipeline/tests/data/output_dir/QCJob/Gerwi"
                         "ck_6123/filtered_sequences/ISB_R1_001.trimmed.fastq."
                         "gz -2 sequence_processing_pipeline/tests/data/output"
                         "_dir/QCJob/Gerwick_6123/filtered_sequences/ISB_R2_00"
                         "1.trimmed.fastq.gz",
            'count': 9},
        'QCJob_NYU_BMS_Melanoma_13059.array-details': {
            'first_line': "fastp --adapter_sequence AACC --adapter_sequence_r2"
                          " GGTT -l 100 -i sequence_processing_pipeline/tests/"
                          "data/output_dir/ConvertJob/NYU_BMS_Melanoma_13059/2"
                          "2_001_710_503_791_00_R1_001.fastq.gz -I sequence_pr"
                          "ocessing_pipeline/tests/data/output_dir/ConvertJob/"
                          "NYU_BMS_Melanoma_13059/22_001_710_503_791_00_R2_001"
                          ".fastq.gz -w 16 -j sequence_processing_pipeline/tes"
                          "ts/data/output_dir/QCJob/NYU_BMS_Melanoma_13059/jso"
                          "n/22_001_710_503_791_00_R1_001.json -h sequence_pro"
                          "cessing_pipeline/tests/data/output_dir/QCJob/NYU_BM"
                          "S_Melanoma_13059/html/22_001_710_503_791_00_R1_001."
                          "html -o sequence_processing_pipeline/tests/data/out"
                          "put_dir/QCJob/NYU_BMS_Melanoma_13059/trimmed_sequen"
                          "ces/22_001_710_503_791_00_R1_001.fastp.fastq.gz -O "
                          "sequence_processing_pipeline/tests/data/output_dir/"
                          "QCJob/NYU_BMS_Melanoma_13059/trimmed_sequences/22_0"
                          "01_710_503_791_00_R2_001.fastp.fastq.gz -R 22_001_7"
                          "10_503_791_00_R1_001_report;fastp --adapter_sequenc"
                          "e AACC --adapter_sequence_r2 GGTT -l 100 -i sequenc"
                          "e_processing_pipeline/tests/data/output_dir/Convert"
                          "Job/NYU_BMS_Melanoma_13059/EP808104A01_R1_001.fastq"
                          ".gz -I sequence_processing_pipeline/tests/data/outp"
                          "ut_dir/ConvertJob/NYU_BMS_Melanoma_13059/EP808104A0"
                          "1_R2_001.fastq.gz -w 16 -j sequence_processing_pipe"
                          "line/tests/data/output_dir/QCJob/NYU_BMS_Melanoma_1"
                          "3059/json/EP808104A01_R1_001.json -h sequence_proce"
                          "ssing_pipeline/tests/data/output_dir/QCJob/NYU_BMS_"
                          "Melanoma_13059/html/EP808104A01_R1_001.html -o sequ"
                          "ence_processing_pipeline/tests/data/output_dir/QCJo"
                          "b/NYU_BMS_Melanoma_13059/trimmed_sequences/EP808104"
                          "A01_R1_001.fastp.fastq.gz -O sequence_processing_pi"
                          "peline/tests/data/output_dir/QCJob/NYU_BMS_Melanoma"
                          "_13059/trimmed_sequences/EP808104A01_R2_001.fastp.f"
                          "astq.gz -R EP808104A01_R1_001_report",
            'last_line': "fastp --adapter_sequence AACC --adapter_sequence_r2 "
                         "GGTT -l 100 -i sequence_processing_pipeline/tests/da"
                         "ta/output_dir/ConvertJob/NYU_BMS_Melanoma_13059/EP80"
                         "5337A01_R1_001.fastq.gz -I sequence_processing_pipel"
                         "ine/tests/data/output_dir/ConvertJob/NYU_BMS_Melanom"
                         "a_13059/EP805337A01_R2_001.fastq.gz -w 16 -j sequenc"
                         "e_processing_pipeline/tests/data/output_dir/QCJob/NY"
                         "U_BMS_Melanoma_13059/json/EP805337A01_R1_001.json -h"
                         " sequence_processing_pipeline/tests/data/output_dir/"
                         "QCJob/NYU_BMS_Melanoma_13059/html/EP805337A01_R1_001"
                         ".html -o sequence_processing_pipeline/tests/data/out"
                         "put_dir/QCJob/NYU_BMS_Melanoma_13059/trimmed_sequenc"
                         "es/EP805337A01_R1_001.fastp.fastq.gz -O sequence_pro"
                         "cessing_pipeline/tests/data/output_dir/QCJob/NYU_BMS"
                         "_Melanoma_13059/trimmed_sequences/EP805337A01_R2_001"
                         ".fastp.fastq.gz -R EP805337A01_R1_001_report",
            'count': 250}
    }


if __name__ == '__main__':
    unittest.main()
