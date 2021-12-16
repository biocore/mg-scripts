from sequence_processing_pipeline.PipelineError import PipelineError
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os.path import exists


class SequenceDirectory:
    def __init__(self, sequence_directory, sample_sheet_path):
        if sequence_directory:
            if exists(sequence_directory):
                self.sequence_directory = sequence_directory
            else:
                raise PipelineError(f"Directory {sequence_directory} "
                                    "does not exist")
        else:
            raise PipelineError("A value for sequence_directory must be "
                                "supplied")

        if sample_sheet_path:
            if exists(sample_sheet_path):
                self.sample_sheet_path = sample_sheet_path
            else:
                raise PipelineError(f"Sample sheet {sample_sheet_path} "
                                    "does not exist")
        else:
            raise PipelineError("An sample sheet must be supplied")

        sheet = KLSampleSheet(sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % sample_sheet_path
            raise PipelineError(s)

        header = valid_sheet.Header
        self.experiment_name = header['Experiment Name']
        self.chemistry = header['Chemistry']
