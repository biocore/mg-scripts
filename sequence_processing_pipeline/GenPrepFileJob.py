from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from os import makedirs, symlink
from os.path import join, exists, basename
from shutil import copytree
from functools import partial
import re


class GenPrepFileJob(Job):
    def __init__(self, run_dir, convert_job_path, qc_job_path, output_path,
                 sample_sheet_path, seqpro_path, project_list, modules_to_load,
                 qiita_job_id, is_amplicon=False):

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
        self.is_amplicon = is_amplicon
        self.prep_file_paths = None

        # make the 'root' of your run_directory
        makedirs(join(self.output_path, self.run_id), exist_ok=True)
        # copy bcl-convert's Stats-equivalent directory to the
        # run_directory
        copytree(join(convert_job_path, 'Reports'),
                 join(self.output_path, self.run_id, 'Reports'))

        for project in project_list:
            src_path = partial(join, qc_job_path, project)
            filtered_seq_dir = src_path('filtered_sequences')
            trimmed_seq_dir = src_path('trimmed_sequences')
            fastp_rept_dir = src_path('fastp_reports_dir', 'json')
            amplicon_seq_dir = join(convert_job_path, project)

            dst = join(self.output_path, self.run_id, project)

            if self.is_amplicon:
                if exists(amplicon_seq_dir):
                    makedirs(dst, exist_ok=True)
                    symlink(amplicon_seq_dir, join(dst, 'amplicon'))
            else:
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
        self.command = [self.seqpro_path, '--verbose',
                        join(self.output_path, self.run_id),
                        f'"{self.sample_sheet_path}"',
                        join(self.output_path, 'PrepFiles')]

    def _get_prep_file_paths(self, stdout):
        tmp_l = stdout.split('\n')
        tmp_l = [x for x in tmp_l if x != '']
        tmp_d = {}

        for line in tmp_l:
            # assume search will always yield a result on legit output.
            qiita_id = re.search(r'\((\d+)\)$', line)[1]

            if qiita_id not in tmp_d:
                tmp_d[qiita_id] = []

            # extract absolute file-path, removing trailing whitespace.
            tmp_d[qiita_id].append(line.replace(f'({qiita_id})', '').strip())

        return tmp_d

    def run(self, callback=None):
        # note that if GenPrepFileJob will be run after QCJob in a Pipeline,
        # and QCJob currently moves its products to the final location. It
        # would be cleaner if it did not do this, but currently that is how
        # it's done. Hence, self.output_directory and the path to run_dir
        # might be different locations than the others.
        results = self._system_call(' '.join(self.command), callback=callback)

        if results['return_code'] != 0:
            raise PipelineError("Seqpro encountered an error")

        # if successful, store results.
        self.prep_file_paths = self._get_prep_file_paths(results['stdout'])
