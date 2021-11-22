from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os import makedirs, symlink
from os.path import join, exists, basename
from shutil import copytree


class GenPrepFileJob(Job):
    def __init__(self, run_dir, convert_job_path, qc_job_path, output_path,
                 sample_sheet_path, seqpro_path, project_list, modules_to_load,
                 qiita_job_id):

        super().__init__(run_dir,
                         output_path,
                         'GenPrepFileJob',
                         [seqpro_path],
                         1000,
                         modules_to_load=modules_to_load)

        self.run_id = basename(run_dir)
        self.sample_sheet_path = sample_sheet_path
        self.seqpro_path = seqpro_path
        self.qiita_job_id = qiita_job_id

        # make the 'root' of your run_directory
        makedirs(join(self.output_path, self.run_id), exist_ok=True)
        # copy bcl-convert's Stats-equivalent directory to the the
        # run_directory
        copytree(join(convert_job_path, 'Reports'),
                 join(self.output_path, self.run_id, 'Reports'))

        for project in project_list:
            filtered_seq_dir = join(qc_job_path, project, 'filtered_sequences')
            trimmed_seq_dir = join(qc_job_path, project, 'trimmed_sequences')
            fastp_rept_dir = join(qc_job_path,
                                  project,
                                  'fastp_reports_dir',
                                  'json')
            dst = join(self.output_path, self.run_id, project)

            if exists(filtered_seq_dir):
                makedirs(dst, exist_ok=True)
                symlink(filtered_seq_dir, join(dst, 'filtered_sequences'))

            if exists(trimmed_seq_dir):
                makedirs(dst, exist_ok=True)
                symlink(trimmed_seq_dir, join(dst, 'trimmed_sequences'))

            if exists(fastp_rept_dir):
                makedirs(dst, exist_ok=True)
                symlink(fastp_rept_dir, join(dst, 'json'))

        # seqpro usage:
        # seqpro path/to/run_dir path/to/sample/sheet /path/to/fresh/output_dir

        # note that seqpro takes in a sample-sheet as input. Unlike other
        # Jobs that process the sample-sheet, validate parameters, and ensure
        # that output is segregated by project name, seqpro will be doing that
        # for us. A single call to seqpro will generate n output files, one
        # for each project described in the sample-sheet's Bioinformatics
        # heading.
        self.command = [self.seqpro_path,
                        join(self.output_path, self.run_id),
                        self.sample_sheet_path,
                        join(self.output_path, 'PrepFiles')]

    def run(self):
        # note that if GenPrepFileJob will be run after QCJob in a Pipeline,
        # and QCJob currently moves its products to the final location. It
        # would be cleaner if it did not do this, but currently that is how
        # it's done. Hence, self.output_directory and the path to run_dir
        # might be different locations than the others.
        results = self._system_call(' '.join(self.command))

        if results['return_code'] != 0:
            raise PipelineError("Seqpro encountered an error")
