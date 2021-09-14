import shutil
from shutil import SameFileError
import os
import logging
from sequence_processing_pipeline.PipelineError import PipelineError
from metapool import KLSampleSheet, validate_sample_sheet
from os.path import join, split, exists
import json


class SequenceDirectory:
    def __init__(self, sequence_directory, external_sample_sheet=None,
                 backup_sample_sheet=False):
        if sequence_directory:
            if exists(sequence_directory):
                self.sequence_directory = sequence_directory
            else:
                s = "Directory %s does not exist." % sequence_directory
                raise PipelineError(s)
        else:
            s = "A value for sequence_directory must be supplied."
            raise PipelineError(s)

        if external_sample_sheet:
            if exists(external_sample_sheet):
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

            # TODO: Disable Fuzzy finder. If the user did not provide a sample sheet,
            #  the user will not get a run out the other end.
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

        s = join(self.sequence_directory, 'Data', 'FastQ')
        self.fastq_results_directory = s

        try:
            os.makedirs(self.fastq_results_directory, mode=750, exist_ok=True)
        except OSError as e:
            # this is a known potential error. Re-raise it as a
            # PinelineError, so it gets handled in the same location as the
            # others.
            raise PipelineError(str(e))

        if backup_sample_sheet:
            self.sheet_backup_dir = join(self.sequence_directory,
                                                 "orig_sample_sheets")
            os.makedirs(self.sheet_backup_dir, mode=750, exist_ok=True)

        # if not exists(join_path(self.seq_dir, 'alockfile')):
        # raise PipelineError("data already processing? %s" % self.seq_dir)
        # if exists(join_path(self.seq_dir, 'processed')):
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

        # extract needed metadata from the sample sheet.
        results = self._process_sample_sheet(input)

        # add to results other valuable metadata
        results['sequence_directory'] = self.sequence_directory

        if self.external_sample_sheet:
            results['sample_sheet_path'] = self.external_sample_sheet
        else:
            results['sample_sheet_path'] = self.sample_sheets[0]

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
                    csv_list.append(join(root, some_file))

        return csv_list

    def _backup_sample_sheet(self, csv_file_path):
        '''
        Store a copy of the selected sample sheet in
        self.sample_sheet_backup_directory.
        :param csv_file_path: path to the sample sheet
        '''

        # copy orig sample sheet to separate directory. This does NOT create
        # numerical increments of file copies. One copy only of latest file.
        file_name = split(csv_file_path)[1] + '.bak'

        archived_path = join(self.sheet_backup_dir, file_name)

        logging.debug("BACKUP DIR: %s" % self.sheet_backup_dir)
        logging.debug("ARCHIVED_PATH: %s" % archived_path)
        logging.debug("SAMPLE SHEET PATH: %s" % csv_file_path)

        try:
            shutil.copyfile(csv_file_path, archived_path)
        except OSError as e:
            # re-raise known possible errors as PipelineErrors.
            raise PipelineError(str(e))
        except SameFileError as e:
            raise PipelineError(str(e))

    def _process_sample_sheet(self, sample_sheet):
        results = {}
        sheet = KLSampleSheet(sample_sheet)
        valid_sheet = validate_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % sample_sheet
            raise PipelineError(s)

        header = valid_sheet.Header
        experiment_name = header['Experiment Name']
        chemistry = header['Chemistry']
        data = header['Data']
        reads = valid_sheet.Reads
        I7_Index_ID = valid_sheet.samples[0]['I7_Index_ID']
        job_index_val = 'I12'

        if len(reads) == 1:
            # direction == 1
            rn1 = reads[0]
            base_mask = "--use-bases-mask Y{},{}".format(rn1, job_index_val)
        elif len(reads) == 2:
            # direction == 2
            rn1 = reads[0]
            base_mask = "--use-bases-mask Y{},{},{},Y{}".format(rn1,
                                                                job_index_val,
                                                                job_index_val,
                                                                rn1)
        else:
            raise PipelineError("Unexpected number of reads: %s" % str(reads))

        results['chemistry'] = chemistry
        results['base_mask'] = base_mask
        results['experiment_name'] = experiment_name

        l = []
        for item in json.loads(valid_sheet.to_json())['Data']:
            d = {}
            if 'Sample_Project' in item:
                d['project_name'] = item['Sample_Project']
            elif 'Project' in item:
                d['project_name'] = item['Project']
            else:
                raise PipelineError("Cannot determine project column in %s" % sample_sheet)

            l.append(d)

        results['data'] = l

        return results