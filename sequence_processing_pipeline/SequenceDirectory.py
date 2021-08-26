import shutil
import pandas as pd
import os
from sequence_processing_pipeline.exceptions import PipelineError
import logging


class SequenceDirectory:
    def __init__(self, seq_dir, rta_file_name):
        logging.debug("Creating SequenceDirectory Object")
        self.seq_dir = seq_dir
        self.check_file = os.path.join(self.seq_dir, rta_file_name)
        self.csv_list = []
        self.mk_path = None
        self.sample_sheet_backup_directory = os.path.join(self.seq_dir, "orig_sample_sheets")
        self.map_n_count = {'12': 'I12', '8': 'I8', '0': 'I12'}

    def _get_csv_list(self):
        '''
        Return a list of CSVs found in the directory.
        :return:
        '''
        for root, dirs, files in os.walk(self.seq_dir):
            for some_file in files:
                # look for files that may be in upper or mixed case as well.
                if some_file.lower().endswith('.csv'):
                    self.csv_list.append(os.path.join(root, some_file))

    def _perform_sanity_checks(self):
        '''

        :return:
        '''
        if not os.path.exists(self.check_file):
            # this error isn't likely to occur as Pipeline object checks for
            # this before creating an SequenceDirectory object.
            raise PipelineError("RTAComplete file doesn't exist.")

        if not self.csv_list:
            # this is a fatal (with respect to this directory) error, as
            # there's no work to do.
            raise PipelineError("RTAComplete file present, but sample sheets are not present for %s." % self.seq_dir)

        # if not os.path.exists(os.path.join(self.seq_dir, 'alockfile')):
            # raise PipelineError("data already processing? %s" % self.seq_dir)

        if os.path.exists(os.path.join(self.seq_dir, 'processed')):
            raise PipelineError("bcl conversion complete? %s" % self.seq_dir)

    def prepare_data_location(self):
        '''
        Look for sample-sheets, perform sanity checks, and make directories.
        :return:
        '''
        try:
            self._get_csv_list()
            self._perform_sanity_checks()
            self.mk_path = os.path.join(self.seq_dir, "Data/Fastq")
            os.makedirs(self.mk_path, mode=750, exist_ok=True)
        except Exception as e:
            logging.error(str(e))
            # TODO: Send email notification
            return False

        # data locations have been properly vetted and directories have been
        # made.
        return True

    def _backup_sample_sheet(self, csv_file_path):
        '''

        :param csv_file_path:
        :return:
        '''
        # we don't want to create this directory on init() if we have to early
        # abort. We should create it at the proper time. In which case, we'll
        # test to see if the path exists each time, before we copy.
        if not os.path.exists(self.sample_sheet_backup_directory):
            # just in case there is more than one level being created, we'll
            # use makedirs().
            # TODO RESTORE mode=750 after testing os.makedirs(self.sample_sheet_backup_directory, mode=750, exist_ok=True)
            os.makedirs(self.sample_sheet_backup_directory, exist_ok=True)

        logging.debug("BACKUP DIR: %s" % self.sample_sheet_backup_directory)
        # copy orig sample sheet to separate directory. This does NOT create
        # numerical increments of file copies. One copy only of latest file.
        file_name_only = os.path.split(csv_file_path)[1]
        archived_path = os.path.join(self.sample_sheet_backup_directory, file_name_only + '.bak')

        logging.debug("ARCHIVED_PATH: %s" % archived_path)
        logging.debug("CSV FILE PATH: %s" % csv_file_path)

        # TODO: Consider handling the exception raised if this fails.
        shutil.copyfile(csv_file_path, archived_path)

    def _process_single_file(self, csv_file):
        '''

        :param csv_file:
        :return:
        '''
        logging.debug("Reading %s" % csv_file)
        df = pd.read_csv(csv_file, header=None)

        # drop all lines with no data (all "NaN"). Compresses sample sheet to
        # no blank lines.
        df = df.dropna(how='all')
        df = df.fillna(".")

        # Personally I feel initializing them to None is a better practice,
        # but we currently we rely on 0 to be the default value if the
        # elements are not found in the sheet.
        contact_index = 0
        bio_index = 0
        info_dict = {}

        # proposed default values for read_index 1 and 2, if values are not
        # found.
        read_index_1 = 151
        read_index_2 = 151

        for index in range(len(df)):
            if 'Bioinformatics' in df.iloc[index, 0]:
                bio_index = index
                print("bio_index=", index)
            elif '[Contact' in df.iloc[index, 0]:
                print("contact_index=", index)
                contact_index = index + 1
            else:
                contact_index = len(df)

            if '[Reads]' in df.iloc[index, 0]:
                read_index_1 = index
                read_index_2 = index + 3

            # assemble dictionary for csv variables

            keyword_list = ['Experiment', 'Assay', 'Chemistry', 'ReverseComplement']
            info_df = pd.DataFrame()

            for keyword in keyword_list:
                if keyword in df.iloc[index, 0]:
                    info_df = info_df.append(df.iloc[index])
                    value = df.iloc[index, 1]
                    info_dict[keyword] = value

        if info_dict["Chemistry"] == "Amplicon":
            logging.debug("Amplicon chemistry true. Removing false barcodes")

            n_count = None
            for line in csv_file:
                for entry in ["NNNNNNNN", "NNNNNNNNNNNN"]:
                    if entry in line:
                        n_count = str(len(entry))
                        logging.debug("N_count: %s" % n_count)
                        # assume once we find it in this file, we don't
                        # need to look further, as that would just
                        # overwrite the previous value.
                        break

            # just in case there are instances where n_count won't be found,
            # this preserves the original idea that 'I12' is the default.
            job_index_val = 'I12'
            if n_count:
                if n_count in self.map_n_count:
                    # if n_count found does not match a known value, we will
                    # default to 'I12'. If we should raise an Error instead,
                    # we should do that.
                    job_index_val = self.map_n_count[n_count]

            base_mask = "--use-bases-mask Y150,%s,%s,Y150" % (job_index_val, job_index_val)

            ### check to see if both Read values are present at 150/151
            ### if both, direction==2
            ### else direction==1
            ### used for bases-mask in bcl2fastq
            ### pull from read_df or develop other method

            ### replace N strings if Amplicon Sample sheets with blank for processing.
            ### can be moved to where we parse chemistry variable.
            ### we have already copied original to backup location.
            ### create new sample sheet with same name as original, removing false barcodes
            csv_1 = open(csv_file, 'r')

            ### should check to see if NNNNNNNN is present but this is historical
            csv_1 = ''.join([i for i in csv_1]).replace("NNNNNNNNNNNN", "")

            # I think instead of making a '.bak' file, reading the original
            # file for processing, and then reopening the file to be
            # overwritten, it would be better to simply read the original
            # file into memory (ala DF), close it, and write out the processed
            # contents of the original file to a new file. The names can
            # always be renamed afterward.
            csv_2 = open(csv_file, 'w')
            csv_2.writelines(csv_1)
            csv_2.close()

        bioinfo_df = pd.DataFrame(df.iloc[bio_index:contact_index - 1])

        logging.debug("contact_index == ", contact_index)

        contact_df = pd.DataFrame(df.iloc[contact_index:])
        if contact_df.empty:
            contact_df = pd.DataFrame({"some_email@gmail.com"})
        logging.debug(contact_df.empty)

        logging.debug("after empty df check")
        read_df = pd.DataFrame(df.iloc[read_index_1:read_index_2])

        logging.debug(bioinfo_df)
        logging.debug(contact_df)

        logging.debug(read_df.T)
        read_df = read_df.dropna(how='all')
        logging.debug(read_df)

    def process_data(self):
        '''

        :return:
        '''
        results = []
        for csv_file in self.csv_list:
            logging.debug("Processing %s..." + csv_file)

            # first back up csv file
            self._backup_sample_sheet(csv_file)

            # then process a single csv file
            # metadata from each processed file is added to a list of results
            # that can then be handed off for submitting jobs.
            results.append(self._process_single_file(csv_file))

        return results