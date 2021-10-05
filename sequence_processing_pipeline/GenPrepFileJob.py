from sequence_processing_pipeline.Job import Job
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
        self.seqpro_path = seqpro_path
        self.modules_to_load = modules_to_load
        self.output_directory = output_directory

    def run(self):
        # seqpro usage:
        # seqpro path/to/run_dir path/to/sample/sheet /path/to/fresh/output_dir
        cmd = [self.seqpro_path]
        cmd.append(self.run_dir)
        cmd.append(self.sample_sheet_path)

        # note that seqpro takes in a sample-sheet as input. Unlike other
        # Jobs that process the sample-sheet, validate parameters, and ensure
        # that output is segregated by project name, seqpro will be doing that
        # for us. Hence our output directory will just be a base directory
        # that seqpro can use to write project directories into.
        tmp = join(self.output_directory, 'FastQC_output')
        cmd.append(tmp)

        # note that if GenPrepFileJob will be run after QCJob in a Pipeline,
        # and QCJob currently moves its products to the final location. It
        # would be cleaner if it did not do this, but currently that is how
        # it's done. Hence, self.output_directory and the path to run_dir
        # might be different locations than the others.

        out, err, rc = self._system_call(' '.join(cmd))

        logging.debug(f"Seqpro STDOUT: {out}")
        logging.debug(f"Seqpro STDERR: {err}")
        logging.debug(f"Seqpro return code: {rc}")

        if rc != 0:
            raise PipelineError("Seqpro encountered an error")
