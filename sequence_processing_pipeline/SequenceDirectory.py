import shutil
from shutil import SameFileError
import pandas as pd
import os
import logging
from sequence_processing_pipeline.PipelineError import PipelineError


class SequenceDirectory:
    def __init__(self, sequence_directory, external_sample_sheet=None,
                 backup_sample_sheet=False):
        if sequence_directory:
            if os.path.exists(sequence_directory):
                self.sequence_directory = sequence_directory
            else:
                s = "Directory %s does not exist." % sequence_directory
                raise PipelineError(s)
        else:
            s = "A value for sequence_directory must be supplied."
            raise PipelineError(s)

        if external_sample_sheet:
            if os.path.exists(external_sample_sheet):
                self.external_sample_sheet = external_sample_sheet
            else:
                s = "External sample sheet {} does not exist."
                raise PipelineError(s.format(self.external_sample_sheet))
        else:
            # if a path to an external sample sheet wasn't provided, SDO will
            # look for sample sheets within the directory itself. Since sample
            # sheets can appear more than once, it will prioritize usage in
            # order of the paths defined in common_locations.
            sample_sheets = self._find_sample_sheets()
            common_locations = ['/SampleSheet.csv',
                                '/Data/Intensities/BaseCalls/SampleSheet.csv']

            prioritized_list = []

            for some_path in sample_sheets:
                relative_path = some_path.replace(self.sequence_directory, '')
                if relative_path in common_locations:
                    logging.info("Sample sheet found: %s" % some_path)
                    prioritized_list.append(some_path)
                else:
                    s = "Sample sheet found in unusual location: {}"
                    logging.warning(s.format(some_path))

            if prioritized_list:
                self.sample_sheets = prioritized_list
            else:
                s = "A suitable sample sheet was not found for %s"
                raise PipelineError(s.format(self.sequence_directory))

        s = os.path.join(self.sequence_directory, 'Data', 'FastQ')
        self.fastq_results_directory = s

        try:
            os.makedirs(self.fastq_results_directory, mode=750, exist_ok=True)
        except OSError as e:
            # this is a known potential error. Re-raise it as a
            # PinelineError, so it gets handled in the same location as the
            # others.
            raise PipelineError(str(e))

        if backup_sample_sheet:
            self.sheet_backup_dir = os.path.join(self.sequence_directory,
                                                 "orig_sample_sheets")
            os.makedirs(self.sheet_backup_dir, mode=750, exist_ok=True)

        # if not os.path.exists(os.path.join(self.seq_dir, 'alockfile')):
        # raise PipelineError("data already processing? %s" % self.seq_dir)
        # if os.path.exists(os.path.join(self.seq_dir, 'processed')):
        # raise PipelineError("bcl conversion complete? %s" % self.seq_dir)

    def process(self):
        if self.external_sample_sheet:
            # if an external sample sheet has been defined, use it as input.
            input = self.external_sample_sheet
        else:
            # use the first sample sheet in the list of sample sheets, as it's
            # in the most correct location. We can assume sample_sheets has
            # been collected already.
            input = self.sample_sheets[0]

        logging.info("Processing %s..." + input)

        # back up sample sheet if asked to.
        if self.sheet_backup_dir:
            self._backup_sample_sheet(input)

        # what this method wants to return to the user is the metadata
        # extracted from the csv file.
        results = self._process_sample_sheet(input)

        results['sequence_directory'] = self.sequence_directory

        if self.external_sample_sheet:
            results['sample_sheet_path'] = self.external_sample_sheet
        else:
            results['sample_sheet_path'] = self.sample_sheets[0]

        results['experiment_name'] = 'GET FROM SAMPLE SHEET OR OTHER SOURCE'

        return results

    def _find_sample_sheets(self):
        '''
        Return a list of sample sheets found in the directory.
        :return:
        '''
        csv_list = []
        for root, dirs, files in os.walk(self.sequence_directory):
            for some_file in files:
                # look for files that may be in upper or mixed case as well.
                if some_file.lower().endswith('.csv'):
                    csv_list.append(os.path.join(root, some_file))

        return csv_list

    def _backup_sample_sheet(self, csv_file_path):
        '''
        Store a copy of the selected sample sheet in
        self.sample_sheet_backup_directory.
        :param csv_file_path: path to the sample sheet
        '''

        # copy orig sample sheet to separate directory. This does NOT create
        # numerical increments of file copies. One copy only of latest file.
        file_name = os.path.split(csv_file_path)[1] + '.bak'

        archived_path = os.path.join(self.sheet_backup_dir, file_name)

        logging.debug("BACKUP DIR: %s" % self.sheet_backup_dir)
        logging.debug("ARCHIVED_PATH: %s" % archived_path)
        logging.debug("CSV FILE PATH: %s" % csv_file_path)

        try:
            shutil.copyfile(csv_file_path, archived_path)
        except OSError as e:
            # re-raise known possible errors as PipelineErrors.
            raise PipelineError(str(e))
        except SameFileError as e:
            raise PipelineError(str(e))

    def _process_sample_sheet(self, csv_file):
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

            keyword_list = ['Experiment', 'Assay', 'Chemistry',
                            'ReverseComplement']
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

            map_n_count = {'12': 'I12', '8': 'I8', '0': 'I12'}

            # just in case there are instances where n_count won't be found,
            # this preserves the original idea that 'I12' is the default.
            job_index_val = 'I12'
            if n_count:
                if n_count in map_n_count:
                    # if n_count found does not match a known value, we will
                    # default to 'I12'. If we should raise an Error instead,
                    # we should do that.
                    job_index_val = map_n_count[n_count]

            base_mask = "--use-bases-mask Y150,{},{},Y150"
            base_mask = base_mask.format(job_index_val, job_index_val)

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
