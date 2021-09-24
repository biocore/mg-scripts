import os
import logging
from sequence_processing_pipeline.PipelineError import PipelineError
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os.path import join, exists


class SequenceDirectory:
    def __init__(self, sequence_directory, external_sample_sheet):
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
                raise PipelineError(s.format(external_sample_sheet))
        else:
            s = "An external sample sheet must be supplied."
            raise PipelineError(s.format(self.external_sample_sheet))

        s = join(self.sequence_directory, 'Data', 'Fastq')
        self.fastq_results_directory = s

        try:
            os.makedirs(self.fastq_results_directory, exist_ok=True)
        except OSError as e:
            # this is a known potential error. Re-raise it as a
            # PipelineError, so it gets handled in the same location as the
            # others.
            raise PipelineError(str(e))

    def process(self):
        logging.info("Processing %s..." % self.external_sample_sheet)

        # extract needed metadata from the sample sheet.
        results = self._process_sample_sheet(self.external_sample_sheet)

        # add to results other valuable metadata
        results['sequence_directory'] = self.sequence_directory
        results['sample_sheet_path'] = self.external_sample_sheet

        return results

    def _process_sample_sheet(self, sample_sheet):
        results = {}
        sheet = KLSampleSheet(sample_sheet)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % sample_sheet
            raise PipelineError(s)

        header = valid_sheet.Header
        experiment_name = header['Experiment Name']
        chemistry = header['Chemistry']

        reads = valid_sheet.Reads
        # I7_Index_ID = valid_sheet.samples[0]['I7_Index_ID']
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

        return results
