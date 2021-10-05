from sequence_processing_pipeline.Job import Job
# from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join
import logging


class GenPrepFileJob(Job):
    def __init__(self, run_dir, sample_sheet_path, output_directory,
                 seqpro_path, modules_to_load):
        # for now, keep this run_dir instead of abspath(run_dir)
        self.job_name = 'GenPrepFileJob'
        super().__init__(run_dir,
                         self.job_name,
                         [seqpro_path],
                         modules_to_load)

        self.sample_sheet_path = sample_sheet_path
        self.products_dir = join(self.run_dir, 'products', 'prep_files')
        self.seqpro_path = seqpro_path
        self.modules_to_load = modules_to_load
        self.destination_directory = output_directory

    def run(self):
        # seqpro path/to/run_dir path/to/sample/sheet /path/to/fresh/output_dir
        cmd = [self.seqpro_path]
        cmd.append(self.run_dir)
        cmd.append(self.sample_sheet_path)
        cmd.append(self.products_dir)

        out, err, rc = self._system_call(' '.join(cmd))

        logging.debug(f"Seqpro STDOUT: {out}")
        logging.debug(f"Seqpro STDERR: {err}")
        logging.debug(f"Seqpro return code: {rc}")

        if rc != 0:
            raise PipelineError("Seqpro encountered an error")
