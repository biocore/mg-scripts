import os


class SlurmBatch:
    def __init__(self, submission_dict):
        self.submission_dict = submission_dict
        pass

    def submit(self):
        '''
        # job submission building
        # mk_path contains fastq output location.
        # assume self.mk_path exists at this point, because we checked and created it earlier.
        # sbatch_submit(mk_path, contact_df, read_df, bioinfo_df, info_dict, csv_file, base_mask, bclconvert_template, dependent_job, experiment_name)
        :return:
        '''
        pass
