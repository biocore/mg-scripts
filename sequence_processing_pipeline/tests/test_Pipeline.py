import json
import os
import shutil

from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
import unittest
from os import utime, makedirs
from os.path import basename, join, abspath
from copy import deepcopy
from time import time
from functools import partial
import re


class TestPipeline(unittest.TestCase):
    def setUp(self):
        package_root = abspath('./sequence_processing_pipeline')
        self.path = partial(join, package_root, 'tests', 'data')
        self.good_config_file = join(package_root, 'configuration.json')
        self.bad_config_file = self.path('bad_configuration.json')
        self.invalid_config_file = 'does/not/exist/configuration.json'
        self.good_run_id = '211021_A00000_0000_SAMPLE'
        self.invalid_run_id = 'not-sample-sequence-directory'
        self.good_output_file_path = self.path('output_dir')
        makedirs(self.good_output_file_path, exist_ok=True)
        self.maxDiff = None
        self.good_sample_sheet_path = self.path('good-sample-sheet.csv')
        self.bad_sample_sheet_path = self.path('duplicate_sample-sample-sheet'
                                               '.csv')
        self.runinfo_file = self.path(self.good_run_id, 'RunInfo.xml')
        self.rtacomplete_file = self.path(self.good_run_id, 'RTAComplete.txt')

        # most of the tests here were written with the assumption that these
        # files already exist.
        self.create_runinfo_file()
        self.create_rtacomplete_file()

        self.sample_sheet_path = ('sequence_processing_pipeline/tests/data/goo'
                                  'd-sample-sheet.csv')
        self.convert_job_working_path = ('/Users/ccowart/Development/mg-script'
                                         's/sequence_processing_pipeline/tests'
                                         '/data/MyConvertJob')
        self.qc_job_working_path = ('/Users/ccowart/Development/mg-scripts/seq'
                                    'uence_processing_pipeline/tests/data/MyQC'
                                    'Job')
        self.fastqc_job_working_path = ('/Users/ccowart/Development/mg-scripts'
                                        '/sequence_processing_pipeline/tests/d'
                                        'ata/MyFastQCJob')

        # create a new instance of the Data/Fastq directory
        base_path = join(self.convert_job_working_path, 'Data', 'Fastq')
        makedirs(join(base_path, 'Feist_11661'), exist_ok=True)
        makedirs(join(base_path, 'Gerwick_6123'), exist_ok=True)
        makedirs(join(base_path, 'NYU_BMS_Melanoma_13059'), exist_ok=True)

        lines = [
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_25_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_58_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_17_64_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_25_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_55_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_19_63_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_24_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_57_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_1_69_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_23_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_46_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_21_51_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_25_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_49_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_29_57_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_24_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_42_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_31_62_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_21_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_41_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_33_50_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_23_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_50_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_3_61_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_22_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_36_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_5_46_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_23_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_41_R2.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51_R1.fastq.gz',
            'AB5075_AZM_TALE_in_MHB_A_baumannii_AB5075_WT_7_51_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-143_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-143_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-144_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-144_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-145_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-145_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-146_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-146_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-147_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-147_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-148_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-148_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-149_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-149_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-150_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-150_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-151_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-151_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-152_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-152_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-153_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-153_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-154_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-154_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-155_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-155_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-156_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-156_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-157_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-157_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-158_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-158_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-159_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-159_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-160_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-160_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-161_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-161_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-162_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-162_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-163_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-163_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-164_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-164_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-165_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-165_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-166_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-166_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-167_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-167_R2.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-168_R1.fastq.gz',
            'CDPH-SAL_Salmonella_Typhi_MDL-168_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_13_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_13_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_28_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_28_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_51_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_10_51_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_14_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_14_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_23_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_23_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_48_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_4_48_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_6_21_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_6_21_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_6_35_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_BOP27_6_35_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_19_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_19_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_35_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_35_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_59_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_18_59_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_16_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_16_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_43_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_43_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_71_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_20_71_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_16_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_16_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_28_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_28_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_52_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_22_52_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_24_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_24_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_52_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_52_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_9_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_24_9_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_27_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_27_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_69_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_69_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_6_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_26_6_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_13_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_13_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_28_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_28_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_53_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_28_53_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_22_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_22_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_60_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_60_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_7_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_30_7_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_20_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_20_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_56_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_56_R2.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_6_R1.fastq.gz',
            'Deoxyribose_PALE_ALE__MG1655_Lib4_32_6_R2.fastq.gz',
            'JBI_KHP_HGL_021_R1.fastq.gz', 'JBI_KHP_HGL_021_R2.fastq.gz',
            'JBI_KHP_HGL_022_R1.fastq.gz', 'JBI_KHP_HGL_022_R2.fastq.gz',
            'JBI_KHP_HGL_023_R1.fastq.gz', 'JBI_KHP_HGL_023_R2.fastq.gz',
            'JBI_KHP_HGL_024_R1.fastq.gz', 'JBI_KHP_HGL_024_R2.fastq.gz',
            'JBI_KHP_HGL_025_R1.fastq.gz', 'JBI_KHP_HGL_025_R2.fastq.gz',
            'JBI_KHP_HGL_026_R1.fastq.gz', 'JBI_KHP_HGL_026_R2.fastq.gz',
            'JBI_KHP_HGL_027_R1.fastq.gz', 'JBI_KHP_HGL_027_R2.fastq.gz',
            'JBI_KHP_HGL_028_Amitesh_soxR_R1.fastq.gz',
            'JBI_KHP_HGL_028_Amitesh_soxR_R2.fastq.gz',
            'JBI_KHP_HGL_029_Amitesh_oxyR_R1.fastq.gz',
            'JBI_KHP_HGL_029_Amitesh_oxyR_R2.fastq.gz',
            'JBI_KHP_HGL_030_Amitesh_soxR_oxyR_R1.fastq.gz',
            'JBI_KHP_HGL_030_Amitesh_soxR_oxyR_R2.fastq.gz',
            'JBI_KHP_HGL_031_Amitesh_rpoS_R1.fastq.gz',
            'JBI_KHP_HGL_031_Amitesh_rpoS_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0326_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0326_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0327_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0327_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0328_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0328_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0329_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0329_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0330_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0330_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0352_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0352_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0353_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0353_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0354_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0354_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0355_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0355_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0356_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0356_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0357_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0357_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0364_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0364_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0366_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0366_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0367_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0367_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0368_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0368_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0369_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0369_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0370_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0370_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0371_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0371_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0372_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0372_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0373_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0373_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0374_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0374_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0375_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0375_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0376_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0376_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0377_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0377_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0378_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0378_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0380_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0380_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0381_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0381_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0382_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0382_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0383_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0383_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0384_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0384_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0385_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0385_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0386_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0386_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0387_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0388_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0388_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0389_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0389_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0390_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0390_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0391_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0391_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0392_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0392_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0393_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0393_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0394_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0394_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0395_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0395_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0396_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0396_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0397_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0397_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0398_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0398_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0399_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0399_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0400_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0400_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0401_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0401_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0402_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0402_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0403_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0403_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0404_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0404_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0405_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0405_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0406_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0406_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0407_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0407_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0408_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0408_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0409_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0409_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0417_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0417_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0418_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0418_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0419_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0419_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0420_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0420_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0421_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0421_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0473_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0473_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0474_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0474_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0483_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0483_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0484_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0484_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0485_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0485_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0486_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0486_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0516_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0516_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0517_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0517_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0518_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0518_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0519_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0519_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0520_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0520_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0521_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0521_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0522_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0522_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0523_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0523_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0524_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0524_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0525_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-B0525_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08624_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08624_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08704_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R08704_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R10727_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11044_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11044_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11078_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11078_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11101_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11101_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11102_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11102_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11103_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11103_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11135_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11135_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11153_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11153_R2.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11154_R1.fastq.gz',
            'JM-MEC__Staphylococcus_aureusstrain_BERTI-R11154_R2.fastq.gz',
            'JM-Metabolic__GN02424_R1.fastq.gz',
            'JM-Metabolic__GN02424_R2.fastq.gz',
            'JM-Metabolic__GN02446_R1.fastq.gz',
            'JM-Metabolic__GN02446_R2.fastq.gz',
            'JM-Metabolic__GN02449_R1.fastq.gz',
            'JM-Metabolic__GN02449_R2.fastq.gz',
            'JM-Metabolic__GN02487_R1.fastq.gz',
            'JM-Metabolic__GN02487_R2.fastq.gz',
            'JM-Metabolic__GN02501_R1.fastq.gz',
            'JM-Metabolic__GN02501_R2.fastq.gz',
            'JM-Metabolic__GN02514_R1.fastq.gz',
            'JM-Metabolic__GN02514_R2.fastq.gz',
            'JM-Metabolic__GN02529_R1.fastq.gz',
            'JM-Metabolic__GN02529_R2.fastq.gz',
            'JM-Metabolic__GN02531_R1.fastq.gz',
            'JM-Metabolic__GN02531_R2.fastq.gz',
            'JM-Metabolic__GN02567_R1.fastq.gz',
            'JM-Metabolic__GN02567_R2.fastq.gz',
            'JM-Metabolic__GN02590_R1.fastq.gz',
            'JM-Metabolic__GN02590_R2.fastq.gz',
            'JM-Metabolic__GN02657_R1.fastq.gz',
            'JM-Metabolic__GN02657_R2.fastq.gz',
            'JM-Metabolic__GN02748_R1.fastq.gz',
            'JM-Metabolic__GN02748_R2.fastq.gz',
            'JM-Metabolic__GN02766_R1.fastq.gz',
            'JM-Metabolic__GN02766_R2.fastq.gz',
            'JM-Metabolic__GN02769_R1.fastq.gz',
            'JM-Metabolic__GN02769_R2.fastq.gz',
            'JM-Metabolic__GN02787_R1.fastq.gz',
            'JM-Metabolic__GN02787_R2.fastq.gz',
            'JM-Metabolic__GN03132_R1.fastq.gz',
            'JM-Metabolic__GN03132_R2.fastq.gz',
            'JM-Metabolic__GN03218_R1.fastq.gz',
            'JM-Metabolic__GN03218_R2.fastq.gz',
            'JM-Metabolic__GN03252_R1.fastq.gz',
            'JM-Metabolic__GN03252_R2.fastq.gz',
            'JM-Metabolic__GN03409_R1.fastq.gz',
            'JM-Metabolic__GN03409_R2.fastq.gz',
            'JM-Metabolic__GN04014_R1.fastq.gz',
            'JM-Metabolic__GN04014_R2.fastq.gz',
            'JM-Metabolic__GN04094_R1.fastq.gz',
            'JM-Metabolic__GN04094_R2.fastq.gz',
            'JM-Metabolic__GN04255_R1.fastq.gz',
            'JM-Metabolic__GN04255_R2.fastq.gz',
            'JM-Metabolic__GN04306_R1.fastq.gz',
            'JM-Metabolic__GN04306_R2.fastq.gz',
            'JM-Metabolic__GN04428_R1.fastq.gz',
            'JM-Metabolic__GN04428_R2.fastq.gz',
            'JM-Metabolic__GN04488_R1.fastq.gz',
            'JM-Metabolic__GN04488_R2.fastq.gz',
            'JM-Metabolic__GN04540_R1.fastq.gz',
            'JM-Metabolic__GN04540_R2.fastq.gz',
            'JM-Metabolic__GN04563_R1.fastq.gz',
            'JM-Metabolic__GN04563_R2.fastq.gz',
            'JM-Metabolic__GN04612_R1.fastq.gz',
            'JM-Metabolic__GN04612_R2.fastq.gz',
            'JM-Metabolic__GN04665_R1.fastq.gz',
            'JM-Metabolic__GN04665_R2.fastq.gz',
            'JM-Metabolic__GN04682_R1.fastq.gz',
            'JM-Metabolic__GN04682_R2.fastq.gz',
            'JM-Metabolic__GN05002_R1.fastq.gz',
            'JM-Metabolic__GN05002_R2.fastq.gz',
            'JM-Metabolic__GN05109_R1.fastq.gz',
            'JM-Metabolic__GN05109_R2.fastq.gz',
            'JM-Metabolic__GN05128_R1.fastq.gz',
            'JM-Metabolic__GN05128_R2.fastq.gz',
            'JM-Metabolic__GN05367_R1.fastq.gz',
            'JM-Metabolic__GN05367_R2.fastq.gz',
            'JM-Metabolic__GN05377_R1.fastq.gz',
            'JM-Metabolic__GN05377_R2.fastq.gz',
            'JM-Metabolic__GN0_2005_R1.fastq.gz',
            'JM-Metabolic__GN0_2005_R2.fastq.gz',
            'JM-Metabolic__GN0_2007_R1.fastq.gz',
            'JM-Metabolic__GN0_2007_R2.fastq.gz',
            'JM-Metabolic__GN0_2009_R1.fastq.gz',
            'JM-Metabolic__GN0_2009_R2.fastq.gz',
            'JM-Metabolic__GN0_2094_R1.fastq.gz',
            'JM-Metabolic__GN0_2094_R2.fastq.gz',
            'JM-Metabolic__GN0_2099_R1.fastq.gz',
            'JM-Metabolic__GN0_2099_R2.fastq.gz',
            'JM-Metabolic__GN0_2148_R1.fastq.gz',
            'JM-Metabolic__GN0_2148_R2.fastq.gz',
            'JM-Metabolic__GN0_2165_R1.fastq.gz',
            'JM-Metabolic__GN0_2165_R2.fastq.gz',
            'JM-Metabolic__GN0_2169_R1.fastq.gz',
            'JM-Metabolic__GN0_2169_R2.fastq.gz',
            'JM-Metabolic__GN0_2172_R1.fastq.gz',
            'JM-Metabolic__GN0_2172_R2.fastq.gz',
            'JM-Metabolic__GN0_2175_R1.fastq.gz',
            'JM-Metabolic__GN0_2175_R2.fastq.gz',
            'JM-Metabolic__GN0_2183_R1.fastq.gz',
            'JM-Metabolic__GN0_2183_R2.fastq.gz',
            'JM-Metabolic__GN0_2215_R1.fastq.gz',
            'JM-Metabolic__GN0_2215_R2.fastq.gz',
            'JM-Metabolic__GN0_2254_R1.fastq.gz',
            'JM-Metabolic__GN0_2254_R2.fastq.gz',
            'JM-Metabolic__GN0_2277_R1.fastq.gz',
            'JM-Metabolic__GN0_2277_R2.fastq.gz',
            'JM-Metabolic__GN0_2290_R1.fastq.gz',
            'JM-Metabolic__GN0_2290_R2.fastq.gz',
            'JM-Metabolic__GN0_2317_R1.fastq.gz',
            'JM-Metabolic__GN0_2317_R2.fastq.gz',
            'JM-Metabolic__GN0_2337_R1.fastq.gz',
            'JM-Metabolic__GN0_2337_R2.fastq.gz',
            'JM-Metabolic__GN0_2354_R1.fastq.gz',
            'JM-Metabolic__GN0_2354_R2.fastq.gz',
            'JM-Metabolic__GN0_2375_R1.fastq.gz',
            'JM-Metabolic__GN0_2375_R2.fastq.gz',
            'JM-Metabolic__GN0_2380_R1.fastq.gz',
            'JM-Metabolic__GN0_2380_R2.fastq.gz',
            'JM-Metabolic__GN0_2393_R1.fastq.gz',
            'JM-Metabolic__GN0_2393_R2.fastq.gz',
            'JM-Metabolic__GN0_2404_R1.fastq.gz',
            'JM-Metabolic__GN0_2404_R2.fastq.gz',
            'P21_E_coli_ELI344_R1.fastq.gz', 'P21_E_coli_ELI344_R2.fastq.gz',
            'P21_E_coli_ELI345_R1.fastq.gz', 'P21_E_coli_ELI345_R2.fastq.gz',
            'P21_E_coli_ELI347_R1.fastq.gz', 'P21_E_coli_ELI347_R2.fastq.gz',
            'P21_E_coli_ELI348_R1.fastq.gz', 'P21_E_coli_ELI348_R2.fastq.gz',
            'P21_E_coli_ELI349_R1.fastq.gz', 'P21_E_coli_ELI349_R2.fastq.gz',
            'P21_E_coli_ELI350_R1.fastq.gz', 'P21_E_coli_ELI350_R2.fastq.gz',
            'P21_E_coli_ELI351_R1.fastq.gz', 'P21_E_coli_ELI351_R2.fastq.gz',
            'P21_E_coli_ELI352_R1.fastq.gz', 'P21_E_coli_ELI352_R2.fastq.gz',
            'P21_E_coli_ELI353_R1.fastq.gz', 'P21_E_coli_ELI353_R2.fastq.gz',
            'P21_E_coli_ELI354_R1.fastq.gz', 'P21_E_coli_ELI354_R2.fastq.gz',
            'P21_E_coli_ELI355_R1.fastq.gz', 'P21_E_coli_ELI355_R2.fastq.gz',
            'P21_E_coli_ELI357_R1.fastq.gz', 'P21_E_coli_ELI357_R2.fastq.gz',
            'P21_E_coli_ELI358_R1.fastq.gz', 'P21_E_coli_ELI358_R2.fastq.gz',
            'P21_E_coli_ELI359_R1.fastq.gz', 'P21_E_coli_ELI359_R2.fastq.gz',
            'P21_E_coli_ELI361_R1.fastq.gz', 'P21_E_coli_ELI361_R2.fastq.gz',
            'P21_E_coli_ELI362_R1.fastq.gz', 'P21_E_coli_ELI362_R2.fastq.gz',
            'P21_E_coli_ELI363_R1.fastq.gz', 'P21_E_coli_ELI363_R2.fastq.gz',
            'P21_E_coli_ELI364_R1.fastq.gz', 'P21_E_coli_ELI364_R2.fastq.gz',
            'P21_E_coli_ELI365_R1.fastq.gz', 'P21_E_coli_ELI365_R2.fastq.gz',
            'P21_E_coli_ELI366_R1.fastq.gz', 'P21_E_coli_ELI366_R2.fastq.gz',
            'P21_E_coli_ELI367_R1.fastq.gz', 'P21_E_coli_ELI367_R2.fastq.gz',
            'P21_E_coli_ELI368_R1.fastq.gz', 'P21_E_coli_ELI368_R2.fastq.gz',
            'P21_E_coli_ELI369_R1.fastq.gz', 'P21_E_coli_ELI369_R2.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_107_BP6_R1.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_107_BP6_R2.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_108_BP7_R1.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_108_BP7_R2.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_109_BP8_R1.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_109_BP8_R2.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_110_M2_R1.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_110_M2_R2.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_111_M5_R1.fastq.gz',
            'Pputida_JBEI__HGL_Pputida_111_M5_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_145_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_145_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_146_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_146_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_147_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_147_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_148_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_148_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_149_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_149_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_150_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_150_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_151_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_151_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_152_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_152_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_153_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_153_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_154_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_154_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_155_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_155_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_156_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_156_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_157_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_157_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_158_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_158_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_159_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_159_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_160_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_160_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_161_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_161_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_162_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_162_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_163_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_163_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_164_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_164_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_165_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_165_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_166_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_166_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_167_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_167_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_168_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_168_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_169_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_169_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_170_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_170_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_171_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_171_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_172_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_172_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_173_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_173_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_174_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_174_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_175_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_175_R2.fastq.gz',
            'Pputida_PALE__HGL_Pputida_176_R1.fastq.gz',
            'Pputida_PALE__HGL_Pputida_176_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_112_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_112_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_113_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_113_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_114_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_114_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_115_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_115_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_116_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_116_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_117_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_117_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_118_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_118_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_119_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_119_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_120_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_120_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_121_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_121_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_122_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_122_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_123_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_123_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_124_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_124_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_125_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_125_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_126_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_126_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_127_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_127_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_128_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_128_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_129_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_129_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_130_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_130_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_131_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_131_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_132_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_132_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_133_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_133_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_134_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_134_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_135_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_135_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_136_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_136_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_137_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_137_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_138_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_138_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_139_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_139_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_140_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_140_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_141_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_141_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_142_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_142_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_143_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_143_R2.fastq.gz',
            'Pputida_TALE__HGL_Pputida_144_R1.fastq.gz',
            'Pputida_TALE__HGL_Pputida_144_R2.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97D_R1.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97D_R2.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97E_R1.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97E_R2.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97L_R1.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97L_R2.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97N_R1.fastq.gz',
            'RMA_KHP_rpoS_Mage_Q97N_R2.fastq.gz',
            'stALE_E_coli_A10_F131_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A10_F131_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A10_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A10_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A10_F43_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A10_F43_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A11_F119_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A11_F119_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A11_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A11_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A11_F43_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A11_F43_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A12_F136_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A12_F136_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A12_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A12_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A12_F43_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A12_F43_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A13_F121_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A13_F121_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A13_F20_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A13_F20_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A13_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A13_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A14_F133_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A14_F133_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A14_F20_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A14_F20_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A14_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A14_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A15_F117_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A15_F117_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A15_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A15_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A15_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A15_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A16_F134_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A16_F134_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A16_F20_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A16_F20_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A16_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A16_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A17_F118_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A17_F118_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A17_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A17_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A18_F130_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A18_F130_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A18_F18_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A18_F18_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A18_F39_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A18_F39_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A1_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A1_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A2_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A2_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A3_F18_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A3_F18_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A3_F40_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A3_F40_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A4_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A4_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A4_F21_I1_R2_R1.fastq.gz',
            'stALE_E_coli_A4_F21_I1_R2_R2.fastq.gz',
            'stALE_E_coli_A4_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A4_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A5_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A5_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A5_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A5_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A6_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A6_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A6_F43_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A6_F43_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A7_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A7_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A7_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A7_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A8_F20_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A8_F20_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A8_F42_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A8_F42_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A9_F21_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A9_F21_I1_R1_R2.fastq.gz',
            'stALE_E_coli_A9_F44_I1_R1_R1.fastq.gz',
            'stALE_E_coli_A9_F44_I1_R1_R2.fastq.gz']
        lines = [join(base_path, 'Feist_11661', x) for x in lines]
        for line in lines:
            with open(line, 'w') as f2:
                f2.write("This is a file.")

        lines = ['3A_R1.fastq.gz', '3A_R2.fastq.gz', '4A_R1.fastq.gz',
                 '4A_R2.fastq.gz', '5B_R1.fastq.gz', '5B_R2.fastq.gz',
                 '6A_R1.fastq.gz', '6A_R2.fastq.gz', '7A_R1.fastq.gz',
                 '7A_R2.fastq.gz', '8A_R1.fastq.gz', '8A_R2.fastq.gz',
                 'GFR_R1.fastq.gz', 'GFR_R2.fastq.gz', 'ISB_R1.fastq.gz',
                 'ISB_R2.fastq.gz']
        lines = [join(base_path, 'Gerwick_6123', x) for x in lines]
        for line in lines:
            with open(line, 'w') as f2:
                f2.write("This is a file.")

        lines = ['22_001_710_503_791_00_R1.fastq.gz',
                 '22_001_710_503_791_00_R2.fastq.gz',
                 '22_001_801_552_503_00_R1.fastq.gz',
                 '22_001_801_552_503_00_R2.fastq.gz',
                 'AP006367B02_R1.fastq.gz',
                 'AP006367B02_R2.fastq.gz', 'AP029018B01_R1.fastq.gz',
                 'AP029018B01_R2.fastq.gz', 'AP032412B04_R1.fastq.gz',
                 'AP032412B04_R2.fastq.gz', 'AP032413B04_R1.fastq.gz',
                 'AP032413B04_R2.fastq.gz', 'AP046324B02_R1.fastq.gz',
                 'AP046324B02_R2.fastq.gz', 'AP046327B02_R1.fastq.gz',
                 'AP046327B02_R2.fastq.gz', 'AP062219B03_R1.fastq.gz',
                 'AP062219B03_R2.fastq.gz', 'AP065292B01_R1.fastq.gz',
                 'AP065292B01_R2.fastq.gz', 'AP103463B01_R1.fastq.gz',
                 'AP103463B01_R2.fastq.gz', 'AP173299B04_R1.fastq.gz',
                 'AP173299B04_R2.fastq.gz', 'AP173301B04_R1.fastq.gz',
                 'AP173301B04_R2.fastq.gz', 'AP173305B04_R1.fastq.gz',
                 'AP173305B04_R2.fastq.gz', 'AP223470B01_R1.fastq.gz',
                 'AP223470B01_R2.fastq.gz', 'AP298002B02_R1.fastq.gz',
                 'AP298002B02_R2.fastq.gz', 'AP309872B03_R1.fastq.gz',
                 'AP309872B03_R2.fastq.gz', 'AP324642B04_R1.fastq.gz',
                 'AP324642B04_R2.fastq.gz', 'AP470859B01_R1.fastq.gz',
                 'AP470859B01_R2.fastq.gz', 'AP481403B02_R1.fastq.gz',
                 'AP481403B02_R2.fastq.gz', 'AP531397B04_R1.fastq.gz',
                 'AP531397B04_R2.fastq.gz', 'AP549678B01_R1.fastq.gz',
                 'AP549678B01_R2.fastq.gz', 'AP549681B02_R1.fastq.gz',
                 'AP549681B02_R2.fastq.gz', 'AP568785B04_R1.fastq.gz',
                 'AP568785B04_R2.fastq.gz', 'AP568787B02_R1.fastq.gz',
                 'AP568787B02_R2.fastq.gz', 'AP581451B02_R1.fastq.gz',
                 'AP581451B02_R2.fastq.gz', 'AP616837B04_R1.fastq.gz',
                 'AP616837B04_R2.fastq.gz', 'AP668628B04_R1.fastq.gz',
                 'AP668628B04_R2.fastq.gz', 'AP668631B04_R1.fastq.gz',
                 'AP668631B04_R2.fastq.gz', 'AP687591B04_R1.fastq.gz',
                 'AP687591B04_R2.fastq.gz', 'AP696363B02_R1.fastq.gz',
                 'AP696363B02_R2.fastq.gz', 'AP732307B04_R1.fastq.gz',
                 'AP732307B04_R2.fastq.gz', 'AP744361A02_R1.fastq.gz',
                 'AP744361A02_R2.fastq.gz', 'AP745799A04_R1.fastq.gz',
                 'AP745799A04_R2.fastq.gz', 'AP771472A04_R1.fastq.gz',
                 'AP771472A04_R2.fastq.gz', 'AP780167B02_R1.fastq.gz',
                 'AP780167B02_R2.fastq.gz', 'AP787247B04_R1.fastq.gz',
                 'AP787247B04_R2.fastq.gz', 'AP795068B04_R1.fastq.gz',
                 'AP795068B04_R2.fastq.gz', 'AP891020A04_R1.fastq.gz',
                 'AP891020A04_R2.fastq.gz', 'AP905750A02_R1.fastq.gz',
                 'AP905750A02_R2.fastq.gz', 'AP911328B01_R1.fastq.gz',
                 'AP911328B01_R2.fastq.gz', 'AP953594A02_R1.fastq.gz',
                 'AP953594A02_R2.fastq.gz', 'AP959450A03_R1.fastq.gz',
                 'AP959450A03_R2.fastq.gz', 'AP967057A04_R1.fastq.gz',
                 'AP967057A04_R2.fastq.gz', 'C14_R1.fastq.gz',
                 'C14_R2.fastq.gz', 'C18_R1.fastq.gz', 'C18_R2.fastq.gz',
                 'C20_R1.fastq.gz', 'C20_R2.fastq.gz', 'C3_R1.fastq.gz',
                 'C3_R2.fastq.gz', 'C5_R1.fastq.gz', 'C5_R2.fastq.gz',
                 'C6_R1.fastq.gz', 'C6_R2.fastq.gz', 'C9_R1.fastq.gz',
                 'C9_R2.fastq.gz', 'EP001624B01_R1.fastq.gz',
                 'EP001624B01_R2.fastq.gz', 'EP001625B01_R1.fastq.gz',
                 'EP001625B01_R2.fastq.gz', 'EP012991B03_R1.fastq.gz',
                 'EP012991B03_R2.fastq.gz', 'EP023801B04_R1.fastq.gz',
                 'EP023801B04_R2.fastq.gz', 'EP023808B02_R1.fastq.gz',
                 'EP023808B02_R2.fastq.gz', 'EP032410B02_R1.fastq.gz',
                 'EP032410B02_R2.fastq.gz', 'EP032412B02_R1.fastq.gz',
                 'EP032412B02_R2.fastq.gz', 'EP043583B01_R1.fastq.gz',
                 'EP043583B01_R2.fastq.gz', 'EP054632B01_R1.fastq.gz',
                 'EP054632B01_R2.fastq.gz', 'EP061002B01_R1.fastq.gz',
                 'EP061002B01_R2.fastq.gz', 'EP073160B01_R1.fastq.gz',
                 'EP073160B01_R2.fastq.gz', 'EP073209B02_R1.fastq.gz',
                 'EP073209B02_R2.fastq.gz', 'EP073216B01_R1.fastq.gz',
                 'EP073216B01_R2.fastq.gz', 'EP087938B02_R1.fastq.gz',
                 'EP087938B02_R2.fastq.gz', 'EP090129B04_R1.fastq.gz',
                 'EP090129B04_R2.fastq.gz', 'EP112567B02_R1.fastq.gz',
                 'EP112567B02_R2.fastq.gz', 'EP121011B01_R1.fastq.gz',
                 'EP121011B01_R2.fastq.gz', 'EP121013B01_R1.fastq.gz',
                 'EP121013B01_R2.fastq.gz', 'EP128904B02_R1.fastq.gz',
                 'EP128904B02_R2.fastq.gz', 'EP128910B01_R1.fastq.gz',
                 'EP128910B01_R2.fastq.gz', 'EP159692B04_R1.fastq.gz',
                 'EP159692B04_R2.fastq.gz', 'EP159695B01_R1.fastq.gz',
                 'EP159695B01_R2.fastq.gz', 'EP163771B01_R1.fastq.gz',
                 'EP163771B01_R2.fastq.gz', 'EP182060B03_R1.fastq.gz',
                 'EP182060B03_R2.fastq.gz', 'EP182065B04_R1.fastq.gz',
                 'EP182065B04_R2.fastq.gz', 'EP182243B02_R1.fastq.gz',
                 'EP182243B02_R2.fastq.gz', 'EP182346B04_R1.fastq.gz',
                 'EP182346B04_R2.fastq.gz', 'EP184255B04_R1.fastq.gz',
                 'EP184255B04_R2.fastq.gz', 'EP190307B01_R1.fastq.gz',
                 'EP190307B01_R2.fastq.gz', 'EP202095B04_R1.fastq.gz',
                 'EP202095B04_R2.fastq.gz', 'EP202452B01_R1.fastq.gz',
                 'EP202452B01_R2.fastq.gz', 'EP207036B01_R1.fastq.gz',
                 'EP207036B01_R2.fastq.gz', 'EP207041B01_R1.fastq.gz',
                 'EP207041B01_R2.fastq.gz', 'EP207042B04_R1.fastq.gz',
                 'EP207042B04_R2.fastq.gz', 'EP212214B01_R1.fastq.gz',
                 'EP212214B01_R2.fastq.gz', 'EP216516B04_R1.fastq.gz',
                 'EP216516B04_R2.fastq.gz', 'EP230245B01_R1.fastq.gz',
                 'EP230245B01_R2.fastq.gz', 'EP238034B01_R1.fastq.gz',
                 'EP238034B01_R2.fastq.gz', 'EP244360B01_R1.fastq.gz',
                 'EP244360B01_R2.fastq.gz', 'EP244366B01_R1.fastq.gz',
                 'EP244366B01_R2.fastq.gz', 'EP256644B01_R1.fastq.gz',
                 'EP256644B01_R2.fastq.gz', 'EP256645B01_R1.fastq.gz',
                 'EP256645B01_R2.fastq.gz', 'EP260543B04_R1.fastq.gz',
                 'EP260543B04_R2.fastq.gz', 'EP260544B04_R1.fastq.gz',
                 'EP260544B04_R2.fastq.gz', 'EP273332B04_R1.fastq.gz',
                 'EP273332B04_R2.fastq.gz', 'EP282107B01_R1.fastq.gz',
                 'EP282107B01_R2.fastq.gz', 'EP282108B01_R1.fastq.gz',
                 'EP282108B01_R2.fastq.gz', 'EP282276B04_R1.fastq.gz',
                 'EP282276B04_R2.fastq.gz', 'EP291979B04_R1.fastq.gz',
                 'EP291979B04_R2.fastq.gz', 'EP291980B04_R1.fastq.gz',
                 'EP291980B04_R2.fastq.gz', 'EP305735B04_R1.fastq.gz',
                 'EP305735B04_R2.fastq.gz', 'EP316863B03_R1.fastq.gz',
                 'EP316863B03_R2.fastq.gz', 'EP320438B01_R1.fastq.gz',
                 'EP320438B01_R2.fastq.gz', 'EP333541B04_R1.fastq.gz',
                 'EP333541B04_R2.fastq.gz', 'EP337325B04_R1.fastq.gz',
                 'EP337325B04_R2.fastq.gz', 'EP337425B01_R1.fastq.gz',
                 'EP337425B01_R2.fastq.gz', 'EP339053B02_R1.fastq.gz',
                 'EP339053B02_R2.fastq.gz', 'EP339057B02_R1.fastq.gz',
                 'EP339057B02_R2.fastq.gz', 'EP339059B02_R1.fastq.gz',
                 'EP339059B02_R2.fastq.gz', 'EP339061B02_R1.fastq.gz',
                 'EP339061B02_R2.fastq.gz', 'EP372981B04_R1.fastq.gz',
                 'EP372981B04_R2.fastq.gz', 'EP379938B01_R1.fastq.gz',
                 'EP379938B01_R2.fastq.gz', 'EP385379B01_R1.fastq.gz',
                 'EP385379B01_R2.fastq.gz', 'EP385384B01_R1.fastq.gz',
                 'EP385384B01_R2.fastq.gz', 'EP385387B01_R1.fastq.gz',
                 'EP385387B01_R2.fastq.gz', 'EP393712B02_R1.fastq.gz',
                 'EP393712B02_R2.fastq.gz', 'EP393714B01_R1.fastq.gz',
                 'EP393714B01_R2.fastq.gz', 'EP393715B01_R1.fastq.gz',
                 'EP393715B01_R2.fastq.gz', 'EP393717B01_R1.fastq.gz',
                 'EP393717B01_R2.fastq.gz', 'EP393718B01_R1.fastq.gz',
                 'EP393718B01_R2.fastq.gz', 'EP400447B04_R1.fastq.gz',
                 'EP400447B04_R2.fastq.gz', 'EP400448B04_R1.fastq.gz',
                 'EP400448B04_R2.fastq.gz', 'EP410041B01_R1.fastq.gz',
                 'EP410041B01_R2.fastq.gz', 'EP410042B01_R1.fastq.gz',
                 'EP410042B01_R2.fastq.gz', 'EP410046B01_R1.fastq.gz',
                 'EP410046B01_R2.fastq.gz', 'EP422407B01_R1.fastq.gz',
                 'EP422407B01_R2.fastq.gz', 'EP431562B04_R1.fastq.gz',
                 'EP431562B04_R2.fastq.gz', 'EP431570B01_R1.fastq.gz',
                 'EP431570B01_R2.fastq.gz', 'EP431575B01_R1.fastq.gz',
                 'EP431575B01_R2.fastq.gz', 'EP446602B01_R1.fastq.gz',
                 'EP446602B01_R2.fastq.gz', 'EP446604B03_R1.fastq.gz',
                 'EP446604B03_R2.fastq.gz', 'EP446610B02_R1.fastq.gz',
                 'EP446610B02_R2.fastq.gz', 'EP447926B04_R1.fastq.gz',
                 'EP447926B04_R2.fastq.gz', 'EP447927B04_R1.fastq.gz',
                 'EP447927B04_R2.fastq.gz', 'EP447928B04_R1.fastq.gz',
                 'EP447928B04_R2.fastq.gz', 'EP447929B04_R1.fastq.gz',
                 'EP447929B04_R2.fastq.gz', 'EP447940B04_R1.fastq.gz',
                 'EP447940B04_R2.fastq.gz', 'EP447975B02_R1.fastq.gz',
                 'EP447975B02_R2.fastq.gz', 'EP448041B04_R1.fastq.gz',
                 'EP448041B04_R2.fastq.gz', 'EP451428B04_R1.fastq.gz',
                 'EP451428B04_R2.fastq.gz', 'EP455757B04_R1.fastq.gz',
                 'EP455757B04_R2.fastq.gz', 'EP455759B04_R1.fastq.gz',
                 'EP455759B04_R2.fastq.gz', 'EP455763B04_R1.fastq.gz',
                 'EP455763B04_R2.fastq.gz', 'EP479266B04_R1.fastq.gz',
                 'EP479266B04_R2.fastq.gz', 'EP479270B03_R1.fastq.gz',
                 'EP479270B03_R2.fastq.gz', 'EP479794B02_R1.fastq.gz',
                 'EP479794B02_R2.fastq.gz', 'EP479894B04_R1.fastq.gz',
                 'EP479894B04_R2.fastq.gz', 'EP483291B04_R1.fastq.gz',
                 'EP483291B04_R2.fastq.gz', 'EP484973B04_R1.fastq.gz',
                 'EP484973B04_R2.fastq.gz', 'EP487995B04_R1.fastq.gz',
                 'EP487995B04_R2.fastq.gz', 'EP504030B04_R1.fastq.gz',
                 'EP504030B04_R2.fastq.gz', 'EP529635B02_R1.fastq.gz',
                 'EP529635B02_R2.fastq.gz', 'EP533388B01_R1.fastq.gz',
                 'EP533388B01_R2.fastq.gz', 'EP533389B03_R1.fastq.gz',
                 'EP533389B03_R2.fastq.gz', 'EP533426B03_R1.fastq.gz',
                 'EP533426B03_R2.fastq.gz', 'EP533429B04_R1.fastq.gz',
                 'EP533429B04_R2.fastq.gz', 'EP542577B04_R1.fastq.gz',
                 'EP542577B04_R2.fastq.gz', 'EP542578B04_R1.fastq.gz',
                 'EP542578B04_R2.fastq.gz', 'EP554501B04_R1.fastq.gz',
                 'EP554501B04_R2.fastq.gz', 'EP554506B04_R1.fastq.gz',
                 'EP554506B04_R2.fastq.gz', 'EP554513B02_R1.fastq.gz',
                 'EP554513B02_R2.fastq.gz', 'EP554515B04_R1.fastq.gz',
                 'EP554515B04_R2.fastq.gz', 'EP554518B04_R1.fastq.gz',
                 'EP554518B04_R2.fastq.gz', 'EP573296B01_R1.fastq.gz',
                 'EP573296B01_R2.fastq.gz', 'EP573310B01_R1.fastq.gz',
                 'EP573310B01_R2.fastq.gz', 'EP573313B01_R1.fastq.gz',
                 'EP573313B01_R2.fastq.gz', 'EP584756B04_R1.fastq.gz',
                 'EP584756B04_R2.fastq.gz', 'EP587475B04_R1.fastq.gz',
                 'EP587475B04_R2.fastq.gz', 'EP587476B04_R1.fastq.gz',
                 'EP587476B04_R2.fastq.gz', 'EP587477B04_R1.fastq.gz',
                 'EP587477B04_R2.fastq.gz', 'EP587478B04_R1.fastq.gz',
                 'EP587478B04_R2.fastq.gz', 'EP606652B04_R1.fastq.gz',
                 'EP606652B04_R2.fastq.gz', 'EP606656B03_R1.fastq.gz',
                 'EP606656B03_R2.fastq.gz', 'EP606662B04_R1.fastq.gz',
                 'EP606662B04_R2.fastq.gz', 'EP606663B04_R1.fastq.gz',
                 'EP606663B04_R2.fastq.gz', 'EP617440B01_R1.fastq.gz',
                 'EP617440B01_R2.fastq.gz', 'EP617441B01_R1.fastq.gz',
                 'EP617441B01_R2.fastq.gz', 'EP617442B01_R1.fastq.gz',
                 'EP617442B01_R2.fastq.gz', 'EP617443B01_R1.fastq.gz',
                 'EP617443B01_R2.fastq.gz', 'EP636802A01_R1.fastq.gz',
                 'EP636802A01_R2.fastq.gz', 'EP649418A02_R1.fastq.gz',
                 'EP649418A02_R2.fastq.gz', 'EP649623A01_R1.fastq.gz',
                 'EP649623A01_R2.fastq.gz', 'EP649653A04_R1.fastq.gz',
                 'EP649653A04_R2.fastq.gz', 'EP649737A03_R1.fastq.gz',
                 'EP649737A03_R2.fastq.gz', 'EP656055A04_R1.fastq.gz',
                 'EP656055A04_R2.fastq.gz', 'EP657260A01_R1.fastq.gz',
                 'EP657260A01_R2.fastq.gz', 'EP657385A04_R1.fastq.gz',
                 'EP657385A04_R2.fastq.gz', 'EP657386A01_R1.fastq.gz',
                 'EP657386A01_R2.fastq.gz', 'EP667743A04_R1.fastq.gz',
                 'EP667743A04_R2.fastq.gz', 'EP675042B01_R1.fastq.gz',
                 'EP675042B01_R2.fastq.gz', 'EP675044A01_R1.fastq.gz',
                 'EP675044A01_R2.fastq.gz', 'EP675075A04_R1.fastq.gz',
                 'EP675075A04_R2.fastq.gz', 'EP683835A01_R1.fastq.gz',
                 'EP683835A01_R2.fastq.gz', 'EP685640B01_R1.fastq.gz',
                 'EP685640B01_R2.fastq.gz', 'EP702221B04_R1.fastq.gz',
                 'EP702221B04_R2.fastq.gz', 'EP718687A04_R1.fastq.gz',
                 'EP718687A04_R2.fastq.gz', 'EP718688A01_R1.fastq.gz',
                 'EP718688A01_R2.fastq.gz', 'EP721390A04_R1.fastq.gz',
                 'EP721390A04_R2.fastq.gz', 'EP724905B01_R1.fastq.gz',
                 'EP724905B01_R2.fastq.gz', 'EP727972A04_R1.fastq.gz',
                 'EP727972A04_R2.fastq.gz', 'EP729433A02_R1.fastq.gz',
                 'EP729433A02_R2.fastq.gz', 'EP729434A01_R1.fastq.gz',
                 'EP729434A01_R2.fastq.gz', 'EP738468A01_R1.fastq.gz',
                 'EP738468A01_R2.fastq.gz', 'EP738469A01_R1.fastq.gz',
                 'EP738469A01_R2.fastq.gz', 'EP749735A07_R1.fastq.gz',
                 'EP749735A07_R2.fastq.gz', 'EP759450A04_R1.fastq.gz',
                 'EP759450A04_R2.fastq.gz', 'EP768164A02_R1.fastq.gz',
                 'EP768164A02_R2.fastq.gz', 'EP768748A04_R1.fastq.gz',
                 'EP768748A04_R2.fastq.gz', 'EP772143A02_R1.fastq.gz',
                 'EP772143A02_R2.fastq.gz', 'EP772145A02_R1.fastq.gz',
                 'EP772145A02_R2.fastq.gz', 'EP784608A01_R1.fastq.gz',
                 'EP784608A01_R2.fastq.gz', 'EP786631A04_R1.fastq.gz',
                 'EP786631A04_R2.fastq.gz', 'EP790019A01_R1.fastq.gz',
                 'EP790019A01_R2.fastq.gz', 'EP790020A02_R1.fastq.gz',
                 'EP790020A02_R2.fastq.gz', 'EP790021A04_R1.fastq.gz',
                 'EP790021A04_R2.fastq.gz', 'EP790023A01_R1.fastq.gz',
                 'EP790023A01_R2.fastq.gz', 'EP805337A01_R1.fastq.gz',
                 'EP805337A01_R2.fastq.gz', 'EP808104A01_R1.fastq.gz',
                 'EP808104A01_R2.fastq.gz', 'EP808105A01_R1.fastq.gz',
                 'EP808105A01_R2.fastq.gz', 'EP808106A01_R1.fastq.gz',
                 'EP808106A01_R2.fastq.gz', 'EP808109A01_R1.fastq.gz',
                 'EP808109A01_R2.fastq.gz', 'EP808110A04_R1.fastq.gz',
                 'EP808110A04_R2.fastq.gz', 'EP808111A03_R1.fastq.gz',
                 'EP808111A03_R2.fastq.gz', 'EP808112A04_R1.fastq.gz',
                 'EP808112A04_R2.fastq.gz', 'EP843906A04_R1.fastq.gz',
                 'EP843906A04_R2.fastq.gz', 'EP846485A01_R1.fastq.gz',
                 'EP846485A01_R2.fastq.gz', 'EP868682A01_R1.fastq.gz',
                 'EP868682A01_R2.fastq.gz', 'EP872341A01_R1.fastq.gz',
                 'EP872341A01_R2.fastq.gz', 'EP876243A04_R1.fastq.gz',
                 'EP876243A04_R2.fastq.gz', 'EP882752A01_R1.fastq.gz',
                 'EP882752A01_R2.fastq.gz', 'EP886422A01_R1.fastq.gz',
                 'EP886422A01_R2.fastq.gz', 'EP890157A02_R1.fastq.gz',
                 'EP890157A02_R2.fastq.gz', 'EP890158A02_R1.fastq.gz',
                 'EP890158A02_R2.fastq.gz', 'EP899038A04_R1.fastq.gz',
                 'EP899038A04_R2.fastq.gz', 'EP905975A04_R1.fastq.gz',
                 'EP905975A04_R2.fastq.gz', 'EP915769A04_R1.fastq.gz',
                 'EP915769A04_R2.fastq.gz', 'EP921593A04_R1.fastq.gz',
                 'EP921593A04_R2.fastq.gz', 'EP921594A04_R1.fastq.gz',
                 'EP921594A04_R2.fastq.gz', 'EP927458A04_R1.fastq.gz',
                 'EP927458A04_R2.fastq.gz', 'EP927459A04_R1.fastq.gz',
                 'EP927459A04_R2.fastq.gz', 'EP927461A04_R1.fastq.gz',
                 'EP927461A04_R2.fastq.gz', 'EP927462A02_R1.fastq.gz',
                 'EP927462A02_R2.fastq.gz', 'EP929277A02_R1.fastq.gz',
                 'EP929277A02_R2.fastq.gz', 'EP940013A01_R1.fastq.gz',
                 'EP940013A01_R2.fastq.gz', 'EP944059A02_R1.fastq.gz',
                 'EP944059A02_R2.fastq.gz', 'EP970001A01_R1.fastq.gz',
                 'EP970001A01_R2.fastq.gz', 'EP970005A01_R1.fastq.gz',
                 'EP970005A01_R2.fastq.gz', 'EP980752B04_R1.fastq.gz',
                 'EP980752B04_R2.fastq.gz', 'EP981129A02_R1.fastq.gz',
                 'EP981129A02_R2.fastq.gz', 'EP987683A01_R1.fastq.gz',
                 'EP987683A01_R2.fastq.gz', 'EP996831B04_R1.fastq.gz',
                 'EP996831B04_R2.fastq.gz', 'LP127767A01_R1.fastq.gz',
                 'LP127767A01_R2.fastq.gz', 'LP127890A01_R1.fastq.gz',
                 'LP127890A01_R2.fastq.gz', 'LP128476A01_R1.fastq.gz',
                 'LP128476A01_R2.fastq.gz', 'LP128479A01_R1.fastq.gz',
                 'LP128479A01_R2.fastq.gz', 'LP128538A01_R1.fastq.gz',
                 'LP128538A01_R2.fastq.gz', 'LP128539A01_R1.fastq.gz',
                 'LP128539A01_R2.fastq.gz', 'LP128540A01_R1.fastq.gz',
                 'LP128540A01_R2.fastq.gz', 'LP128541A01_R1.fastq.gz',
                 'LP128541A01_R2.fastq.gz', 'LP128543A01_R1.fastq.gz',
                 'LP128543A01_R2.fastq.gz', 'LP154981A01_R1.fastq.gz',
                 'LP154981A01_R2.fastq.gz', 'LP154986A01_R1.fastq.gz',
                 'LP154986A01_R2.fastq.gz', 'LP166715A01_R1.fastq.gz',
                 'LP166715A01_R2.fastq.gz', 'LP169879A01_R1.fastq.gz',
                 'LP169879A01_R2.fastq.gz', 'LP191039A01_R1.fastq.gz',
                 'LP191039A01_R2.fastq.gz', 'LP196272A01_R1.fastq.gz',
                 'LP196272A01_R2.fastq.gz', 'SP205732A02_R1.fastq.gz',
                 'SP205732A02_R2.fastq.gz', 'SP205754A01_R1.fastq.gz',
                 'SP205754A01_R2.fastq.gz', 'SP229387A04_R1.fastq.gz',
                 'SP229387A04_R2.fastq.gz', 'SP230380A02_R1.fastq.gz',
                 'SP230380A02_R2.fastq.gz', 'SP230381A01_R1.fastq.gz',
                 'SP230381A01_R2.fastq.gz', 'SP230382A04_R1.fastq.gz',
                 'SP230382A04_R2.fastq.gz', 'SP231628A02_R1.fastq.gz',
                 'SP231628A02_R2.fastq.gz', 'SP231629A02_R1.fastq.gz',
                 'SP231629A02_R2.fastq.gz', 'SP231630A02_R1.fastq.gz',
                 'SP231630A02_R2.fastq.gz', 'SP231631A02_R1.fastq.gz',
                 'SP231631A02_R2.fastq.gz', 'SP232077A04_R1.fastq.gz',
                 'SP232077A04_R2.fastq.gz', 'SP232079A01_R1.fastq.gz',
                 'SP232079A01_R2.fastq.gz', 'SP232114A04_R1.fastq.gz',
                 'SP232114A04_R2.fastq.gz', 'SP232270A02_R1.fastq.gz',
                 'SP232270A02_R2.fastq.gz', 'SP232309A01_R1.fastq.gz',
                 'SP232309A01_R2.fastq.gz', 'SP232310A04_R1.fastq.gz',
                 'SP232310A04_R2.fastq.gz', 'SP232311A04_R1.fastq.gz',
                 'SP232311A04_R2.fastq.gz', 'SP235186A04_R1.fastq.gz',
                 'SP235186A04_R2.fastq.gz', 'SP235189A01_R1.fastq.gz',
                 'SP235189A01_R2.fastq.gz', 'SP246941A01_R1.fastq.gz',
                 'SP246941A01_R2.fastq.gz', 'SP247340A04_R1.fastq.gz',
                 'SP247340A04_R2.fastq.gz', 'SP257517A04_R1.fastq.gz',
                 'SP257517A04_R2.fastq.gz', 'SP257519A04_R1.fastq.gz',
                 'SP257519A04_R2.fastq.gz', 'SP280481A02_R1.fastq.gz',
                 'SP280481A02_R2.fastq.gz', 'SP284095A03_R1.fastq.gz',
                 'SP284095A03_R2.fastq.gz', 'SP284096A02_R1.fastq.gz',
                 'SP284096A02_R2.fastq.gz', 'SP317293A02_R1.fastq.gz',
                 'SP317293A02_R2.fastq.gz', 'SP317297A02_R1.fastq.gz',
                 'SP317297A02_R2.fastq.gz', 'SP331134A04_R1.fastq.gz',
                 'SP331134A04_R2.fastq.gz', 'SP335002A04_R1.fastq.gz',
                 'SP335002A04_R2.fastq.gz', 'SP353893A02_R1.fastq.gz',
                 'SP353893A02_R2.fastq.gz', 'SP365864A04_R1.fastq.gz',
                 'SP365864A04_R2.fastq.gz', 'SP388683A02_R1.fastq.gz',
                 'SP388683A02_R2.fastq.gz', 'SP399724A04_R1.fastq.gz',
                 'SP399724A04_R2.fastq.gz', 'SP404403A02_R1.fastq.gz',
                 'SP404403A02_R2.fastq.gz', 'SP404405A02_R1.fastq.gz',
                 'SP404405A02_R2.fastq.gz', 'SP404409A02_R1.fastq.gz',
                 'SP404409A02_R2.fastq.gz', 'SP404412A02_R1.fastq.gz',
                 'SP404412A02_R2.fastq.gz', 'SP408629A01_R1.fastq.gz',
                 'SP408629A01_R2.fastq.gz', 'SP410793A01_R1.fastq.gz',
                 'SP410793A01_R2.fastq.gz', 'SP410796A02_R1.fastq.gz',
                 'SP410796A02_R2.fastq.gz', 'SP415021A02_R1.fastq.gz',
                 'SP415021A02_R2.fastq.gz', 'SP415023A02_R1.fastq.gz',
                 'SP415023A02_R2.fastq.gz', 'SP415025A01_R1.fastq.gz',
                 'SP415025A01_R2.fastq.gz', 'SP415030A01_R1.fastq.gz',
                 'SP415030A01_R2.fastq.gz', 'SP416130A04_R1.fastq.gz',
                 'SP416130A04_R2.fastq.gz', 'SP453872A01_R1.fastq.gz',
                 'SP453872A01_R2.fastq.gz', 'SP464350A04_R1.fastq.gz',
                 'SP464350A04_R2.fastq.gz', 'SP464352A03_R1.fastq.gz',
                 'SP464352A03_R2.fastq.gz', 'SP471496A04_R1.fastq.gz',
                 'SP471496A04_R2.fastq.gz', 'SP478193A02_R1.fastq.gz',
                 'SP478193A02_R2.fastq.gz', 'SP490298A02_R1.fastq.gz',
                 'SP490298A02_R2.fastq.gz', 'SP491897A02_R1.fastq.gz',
                 'SP491897A02_R2.fastq.gz', 'SP491898A02_R1.fastq.gz',
                 'SP491898A02_R2.fastq.gz', 'SP491900A02_R1.fastq.gz',
                 'SP491900A02_R2.fastq.gz', 'SP491907A02_R1.fastq.gz',
                 'SP491907A02_R2.fastq.gz', 'SP503615A02_R1.fastq.gz',
                 'SP503615A02_R2.fastq.gz', 'SP506933A04_R1.fastq.gz',
                 'SP506933A04_R2.fastq.gz', 'SP511289A02_R1.fastq.gz',
                 'SP511289A02_R2.fastq.gz', 'SP511294A04_R1.fastq.gz',
                 'SP511294A04_R2.fastq.gz', 'SP515443A04_R1.fastq.gz',
                 'SP515443A04_R2.fastq.gz', 'SP515763A04_R1.fastq.gz',
                 'SP515763A04_R2.fastq.gz', 'SP531696A04_R1.fastq.gz',
                 'SP531696A04_R2.fastq.gz', 'SP561451A04_R1.fastq.gz',
                 'SP561451A04_R2.fastq.gz', 'SP573823A04_R1.fastq.gz',
                 'SP573823A04_R2.fastq.gz', 'SP573824A04_R1.fastq.gz',
                 'SP573824A04_R2.fastq.gz', 'SP573843A04_R1.fastq.gz',
                 'SP573843A04_R2.fastq.gz', 'SP573849A04_R1.fastq.gz',
                 'SP573849A04_R2.fastq.gz', 'SP573859A04_R1.fastq.gz',
                 'SP573859A04_R2.fastq.gz', 'SP573860A01_R1.fastq.gz',
                 'SP573860A01_R2.fastq.gz', 'SP577399A02_R1.fastq.gz',
                 'SP577399A02_R2.fastq.gz', 'SP584547A02_R1.fastq.gz',
                 'SP584547A02_R2.fastq.gz', 'SP584551A08_R1.fastq.gz',
                 'SP584551A08_R2.fastq.gz', 'SP612495A04_R1.fastq.gz',
                 'SP612495A04_R2.fastq.gz', 'SP612496A01_R1.fastq.gz',
                 'SP612496A01_R2.fastq.gz', 'SP631994A04_R1.fastq.gz',
                 'SP631994A04_R2.fastq.gz', 'SP640978A02_R1.fastq.gz',
                 'SP640978A02_R2.fastq.gz', 'SP641029A02_R1.fastq.gz',
                 'SP641029A02_R2.fastq.gz', 'SP645141A03_R1.fastq.gz',
                 'SP645141A03_R2.fastq.gz', 'SP681591A04_R1.fastq.gz',
                 'SP681591A04_R2.fastq.gz', 'SP683466A02_R1.fastq.gz',
                 'SP683466A02_R2.fastq.gz', 'SP704319A04_R1.fastq.gz',
                 'SP704319A04_R2.fastq.gz', 'SP754514A04_R1.fastq.gz',
                 'SP754514A04_R2.fastq.gz', 'ep256643b01_R1.fastq.gz',
                 'ep256643b01_R2.fastq.gz', 'lp127896a01_R1.fastq.gz',
                 'lp127896a01_R2.fastq.gz']
        lines = [join(base_path, 'NYU_BMS_Melanoma_13059', x) for x in lines]
        for line in lines:
            with open(line, 'w') as f2:
                f2.write("This is a file.")

    def tearDown(self):
        # Pipeline is now the only class aware of these files, hence they
        # can be deleted at the end of testing.
        self.delete_runinfo_file()
        self.delete_rtacomplete_file()
        shutil.rmtree(self.convert_job_working_path)

    def make_runinfo_file_unreadable(self):
        os.chmod(self.runinfo_file, 0o000)

    def make_runinfo_file_readable(self):
        os.chmod(self.runinfo_file, 0o777)

    def create_runinfo_file(self):
        with open(self.runinfo_file, 'w') as f:
            f.write("")

    def delete_runinfo_file(self):
        try:
            os.remove(self.runinfo_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def create_rtacomplete_file(self):
        with open(self.rtacomplete_file, 'w') as f:
            f.write("")

    def delete_rtacomplete_file(self):
        try:
            os.remove(self.rtacomplete_file)
        except FileNotFoundError:
            # make method idempotent
            pass

    def test_required_file_checks(self):
        # begin this test by deleting the RunInfo.txt file and verifying that
        # Pipeline object will raise an Error.
        self.delete_runinfo_file()
        with self.assertRaisesRegex(PipelineError, "required file 'RunInfo.xml"
                                                   "' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_output_file_path, 'my_qiita_id', None)

        # delete RTAComplete.txt and recreate RunInfo.txt file to verify that
        # an Error is raised when only RTAComplete.txt is missing.
        self.delete_rtacomplete_file()
        self.create_runinfo_file()
        with self.assertRaisesRegex(PipelineError, "required file 'RTAComplete"
                                                   ".txt' is not present."):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_output_file_path, 'my_qiita_id', None)

        # make RunInfo.xml file unreadable and verify that Pipeline object
        # raises the expected Error.
        self.create_rtacomplete_file()
        self.make_runinfo_file_unreadable()
        with self.assertRaisesRegex(PipelineError, "RunInfo.xml is present, bu"
                                                   "t not readable"):
            Pipeline(self.good_config_file, self.good_run_id,
                     self.good_output_file_path, 'my_qiita_id', None)
        self.make_runinfo_file_readable()

    def test_creation(self):
        # Pipeline should assert due to config_file
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.bad_config_file,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        msg = re.sub(r'not a key in .*?/sequence_processing_pipeline',
                     r'not a key in sequence_processing_pipeline',
                     str(e.exception))
        self.assertEqual(msg, "'search_paths' is not a key in "
                              "sequence_processing_pipeline/tests"
                              "/data/bad_configuration.json")

        # Pipeline should assert due to an invalid config file path.
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.invalid_config_file,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        self.assertEqual(str(e.exception), 'does/not/exist/configuration.json '
                                           'does not exist.')

        # Pipeline should assert on config_file = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(None,
                     self.good_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        self.assertEqual(str(e.exception), 'configuration_file_path cannot be '
                                           'None')

        # Pipeline should assert due to invalid_run_id
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     self.invalid_run_id,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        self.assertEqual(str(e.exception), "A run-dir for 'not-sample-sequence"
                                           "-directory' could not be found")

        # Pipeline should assert on run_id = None
        with self.assertRaises(PipelineError) as e:
            Pipeline(self.good_config_file,
                     None,
                     self.good_output_file_path,
                     'my_qiita_id',
                     None)

        with open(join('sequence_processing_pipeline', 'configuration.json'),
                  'r') as f:
            cfg = json.load(f)
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = -1
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_output_file_path,
                         'my_qiita_id',
                         cpy_cfg)

            self.assertEqual(str(e.exception), 'older_than and younger_than '
                                               'cannot be less than zero.')
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['older_than'] = -1
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_output_file_path,
                         'my_qiita_id',
                         cpy_cfg)
            self.assertEqual(str(e.exception), 'older_than and younger_than '
                                               'cannot be less than zero.')
            with self.assertRaises(PipelineError) as e:
                cpy_cfg = deepcopy(cfg)
                cpy_cfg['configuration']['pipeline']['younger_than'] = 20
                cpy_cfg['configuration']['pipeline']['older_than'] = 30
                Pipeline(self.good_config_file,
                         self.good_run_id,
                         self.good_output_file_path,
                         'my_qiita_id',
                         cpy_cfg)
            self.assertEqual(str(e.exception), 'older_than cannot be equal to '
                                               'or less than younger_than.')

    def test_filter_directories_for_time(self):
        # get a good configuration dict from file.
        with open(join('sequence_processing_pipeline', 'configuration.json'),
                  'r') as f:
            cfg = json.load(f)
            # set a range that 211021_A00000_0000_SAMPLE's timestamps must
            # fall within: between 1 and 2 hours.
            cfg['configuration']['pipeline']['younger_than'] = 2
            cfg['configuration']['pipeline']['older_than'] = 1

            # Set directory's timestamp to between one and two hours and
            # verify it is returned by filter_directories_for_time().

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value older than 1 hour ago + 5 min.
            older_than = current_time - (3600 + (5 * 60))
            tp = self.path('211021_A00000_0000_SAMPLE')
            utime(tp, (older_than, older_than))

            pipeline = Pipeline(self.good_config_file, self.good_run_id,
                                self.good_output_file_path, 'my_qiita_id',
                                cfg)

            obs = pipeline.is_within_time_range(tp)
            self.assertEqual(obs, True)

            # Set directory's timestamp to just older than two hours and
            # verify it is not returned by filter_directories_for_time().

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value older than 2 hour ago + 5 min.
            older_than = current_time - (7200 + (5 * 60))
            utime(tp, (older_than, older_than))

            obs = pipeline.is_within_time_range(tp)
            self.assertEqual(obs, False)

            # Set directory's timestamp to just under one hour and
            # verify it is not returned by filter_directories_for_time().

            # get the current time in seconds since the epoch.
            current_time = time()
            # create an epoch time value younger than 1 hour ago.
            older_than = current_time - 3300
            utime(self.path('211021_A00000_0000_SAMPLE'), (older_than,
                                                           older_than))

            obs = pipeline.is_within_time_range(tp)
            self.assertEqual(obs, False)

    def test_sample_sheet_validation(self):
        # test successful validation of a good sample-sheet
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_output_file_path, 'my_qiita_id',
                            None)

        msgs, val_sheet = pipeline.validate(self.good_sample_sheet_path)

        msgs = [str(x) for x in msgs]

        # a successful validation should return an empty list of error
        # messages and a sheet object.
        self.assertEqual(msgs, [])
        self.assertIsNotNone(val_sheet)

        # test unsuccessful validation of a bad sample-sheet
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_output_file_path, 'my_qiita_id',
                            None)

        msgs, val_sheet = pipeline.validate(self.bad_sample_sheet_path)

        # an unsuccessful validation should return a list of one or more
        # metapool ErrorMessages and a value of None for val_sheet.
        self.assertIsNotNone(msgs)
        self.assertEqual(str(msgs[0]), 'ErrorMessage: A sample already exists '
                                       'with lane 1 and sample-id EP479894B04')

    def test_generate_sample_information_files(self):
        # test sample-information-file generation.
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_output_file_path, 'my_qiita_id',
                            None)

        # assume a validated sample-sheet. Also,
        # generate_sample_information_files calls validate() itself to ensure
        # proper operation.
        paths = partial(pipeline.generate_sample_information_files)
        paths = paths(self.good_sample_sheet_path)

        # confirm files exist in the expected location and with the expected
        # filenames.
        obs = [x.split('sequence_processing_pipeline/')[1] for x in paths]
        exp = ['tests/data/output_dir/NYU_BMS_Melanoma_13059_blanks.tsv',
               'tests/data/output_dir/Feist_11661_blanks.tsv',
               'tests/data/output_dir/Gerwick_6123_blanks.tsv']

        # sort the lists to ensure both are in a fixed order.
        obs.sort()
        exp.sort()

        self.assertEqual(obs, exp)

        # confirm files contain the expected number of lines.
        # This is going to be based on the number of samples named 'BLANK*'
        # in good-sample-sheet.csv.
        exp_lines = {'NYU_BMS_Melanoma_13059_blanks.tsv': 33,
                     'Feist_11661_blanks.tsv': 8,
                     'Gerwick_6123_blanks.tsv': 2}

        exp_first_lines = {
            'NYU_BMS_Melanoma_13059_blanks.tsv': 'BLANK1.1A\t2021-10-21\t193\t'
                                                 'Control\tNegative\tSterile w'
                                                 'ater blank\turban biome\tres'
                                                 'earch facility\tsterile wate'
                                                 'r\tmisc environment\tUSA:CA:'
                                                 'San Diego\tBLANK1.1A\t32.5\t'
                                                 '-117.25\tcontrol blank\tmeta'
                                                 'genome\t256318\tBLANK1.1A\tN'
                                                 'YU_BMS_Melanoma\tTRUE\t'
                                                 'UCSD\tFALSE',
            'Feist_11661_blanks.tsv': 'BLANK.40.12G\t2021-10-21\t193\tControl'
                                      '\tNegative\tSterile water blank\turban '
                                      'biome\tresearch facility\tsterile water'
                                      '\tmisc environment\tUSA:CA:San Diego\tB'
                                      'LANK.40.12G\t32.5\t-117.25\tcontrol bla'
                                      'nk\tmetagenome\t256318\tBLANK.40.12G\t'
                                      'Feist\tTRUE\tUCSD\tFALSE',
            'Gerwick_6123_blanks.tsv': 'BLANK.41.12G\t2021-10-21\t193\tControl'
                                       '\tNegative\tSterile water blank\turban'
                                       ' biome\tresearch facility\tsterile wat'
                                       'er\tmisc environment\tUSA:CA:San Diego'
                                       '\tBLANK.41.12G\t32.5\t-117.25\tcontrol'
                                       ' blank\tmetagenome\t256318\tBLANK.41.1'
                                       '2G\tGerwick\tTRUE\tUCSD\tFALSE'
        }

        exp_last_lines = {
            'NYU_BMS_Melanoma_13059_blanks.tsv': 'BLANK4.4H\t2021-10-21\t193\t'
                                                 'Control\tNegative\tSterile w'
                                                 'ater blank\turban biome\tres'
                                                 'earch facility\tsterile wate'
                                                 'r\tmisc environment\tUSA:CA:'
                                                 'San Diego\tBLANK4.4H\t32.5\t'
                                                 '-117.25\tcontrol blank\tmeta'
                                                 'genome\t256318\tBLANK4.4H\tN'
                                                 'YU_BMS_Melanoma\tTRUE\t'
                                                 'UCSD\tFALSE',
            'Feist_11661_blanks.tsv': 'BLANK.43.12H\t2021-10-21\t193\tControl'
                                      '\tNegative\tSterile water blank\turban'
                                      ' biome\tresearch facility\tsterile wat'
                                      'er\tmisc environment\tUSA:CA:San Diego'
                                      '\tBLANK.43.12H\t32.5\t-117.25\tcontrol'
                                      ' blank\tmetagenome\t256318\tBLANK.43.1'
                                      '2H\tFeist\tTRUE\tUCSD\tFALSE',
            'Gerwick_6123_blanks.tsv': 'BLANK.41.12G\t2021-10-21\t193\tContro'
                                       'l\tNegative\tSterile water blank\turb'
                                       'an biome\tresearch facility\tsterile '
                                       'water\tmisc environment\tUSA:CA:San D'
                                       'iego\tBLANK.41.12G\t32.5\t-117.25\tco'
                                       'ntrol blank\tmetagenome\t256318\tBLAN'
                                       'K.41.12G\tGerwick\tTRUE\tUCSD\t'
                                       'FALSE'
        }

        for some_path in paths:
            some_name = basename(some_path)
            with open(some_path, 'r') as f:
                obs_lines = f.readlines()
                self.assertEqual(len(obs_lines), exp_lines[some_name])
                # confirm that each file contains the expected header.
                header = obs_lines[0].strip()
                self.assertEqual(header, '\t'.join(Pipeline.sif_header))
                # confirm that the first line of each file is as expected.
                obs = obs_lines[1].strip()
                exp = exp_first_lines[some_name]
                self.assertEqual(obs, exp)
                # confirm that the last line of each file is as expected.
                obs = obs_lines[-1].strip()
                exp = exp_last_lines[some_name]
                self.assertEqual(obs, exp)

    def test_audit(self):
        pipeline = Pipeline(self.good_config_file, self.good_run_id,
                            self.good_output_file_path, 'my_qiita_id',
                            None)

        # test the proper extraction of metadata from the sample-sheet.
        results = pipeline.extract_metadata(self.sample_sheet_path)

        sample_ids = results['unique-values']['sample-ids']

        found = pipeline.audit_job(sample_ids, self.convert_job_working_path,
                                   'fastq.gz', 'R')

        not_found = list(set(found) ^ set(sample_ids))

        found.sort()
        not_found.sort()

        self.assertEqual(found[0], '22_001_710_503_791_00')
        self.assertEqual(found[-1], 'stALE_E_coli_A9_F44_I1_R1')
        self.assertEqual(743, len(found))

        # All samples should have been found.
        self.assertEqual(0, len(not_found))

        # verify that the entire list of sample-ids found match what's
        # expected. Since the list of expected ids is very small, we'll
        # perform an exact comparison.
        found = pipeline.audit_job(sample_ids,
                                   self.qc_job_working_path,
                                   'fastq.gz',
                                   'S')

        not_found = list(set(found) ^ set(sample_ids))

        found.sort()
        not_found.sort()

        exp = ["3A", "CDPH-SAL_Salmonella_Typhi_MDL-143"]
        self.assertEqual(found, exp)

        # verify that the first and last sample-ids in the sorted list of
        # sample-ids not found are correct, as well as the number of sample-
        # ids present in the list (or not present in the results).
        not_found = list(set(found) ^ set(sample_ids))
        not_found.sort()
        self.assertEqual(not_found[0], '22_001_710_503_791_00')
        self.assertEqual(not_found[-1], 'stALE_E_coli_A9_F44_I1_R1')
        self.assertEqual(741, len(not_found))

        # assume for now that a corresponding zip file exists for each html
        # file found. Assume for now that all html files will be found in a
        # 'filtered_sequences' or 'trimmed_sequences' subdirectory.
        #
        # verify that the entire list of sample-ids found match what's
        # expected. Since the list of expected ids is very small, we'll
        # perform an exact comparison.
        found = pipeline.audit_job(sample_ids, self.fastqc_job_working_path,
                                   '_fastqc.html', 'S')
        found.sort()
        exp = ["3A", "CDPH-SAL_Salmonella_Typhi_MDL-143"]
        self.assertEqual(found, exp)

        # verify that the first and last sample-ids in the sorted list of
        # sample-ids not found are correct, as well as the number of sample-
        # ids present in the list (or not present in the results).
        not_found = list(set(found) ^ set(sample_ids))
        not_found.sort()
        self.assertEqual(not_found[0], '22_001_710_503_791_00')
        self.assertEqual(not_found[-1], 'stALE_E_coli_A9_F44_I1_R1')
        self.assertEqual(741, len(not_found))


if __name__ == '__main__':
    unittest.main()
