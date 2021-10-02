from sequence_processing_pipeline.Job import Job
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from sequence_processing_pipeline.PipelineError import PipelineError
import logging


class GenPrepFileJob(Job):
    def __init__(self, run_dir, sample_sheet_path, output_directory):
        self.job_name = 'GenPrepFileJob'
        super().__init__(run_dir, self.job_name)
        # metadata = self._process_sample_sheet(sample_sheet_path)
        self.sample_sheet_path = sample_sheet_path
        self.output_directory = output_directory

    def run(self):
        logging.debug('generating prep-files w/seqpro')
        # the run directory, project name, and lane for the prep.
        # from prep.py 424: name = run_id + '.' + project + '.' + lane
        self._system_call(f'seqpro {self.sample_sheet_path}'
                          f'{self.output_directory}')

    def _process_sample_sheet(self, sample_sheet_path):
        sheet = KLSampleSheet(sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % sample_sheet_path
            raise PipelineError(s)

        header = valid_sheet.Header
        chemistry = header['chemistry']
        needs_adapter_trimming = ('TRUE' if
                                  header['Assay'] == 'Metagenomics'
                                  else 'FALSE')

        bioinformatics = valid_sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.

        lst = []
        for i in range(0, len(bioinformatics)):
            lst.append({})

        for key in bioinformatics:
            my_series = bioinformatics[key]
            for index, value in my_series.items():
                # convert all manner of positive/negative implying strings
                # into a true boolean value.
                if value.strip().lower() in ['true', 'yes']:
                    value = True
                elif value.strip().lower() in ['false', 'no']:
                    value = False
                lst[index][key] = value

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry, 'projects': lst,
                'needs_adapter_trimming': needs_adapter_trimming}
