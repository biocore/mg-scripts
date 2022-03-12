from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os.path import basename, join
from os import walk
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
import logging
import re


class ConvertJob(Job):
    def __init__(self, run_dir, output_path, sample_sheet_path, queue_name,
                 node_count, nprocs, wall_time_limit, pmem, bcl_tool_path,
                 modules_to_load, qiita_job_id):
        '''
        ConvertJob provides a convenient way to run bcl-convert or bcl2fastq
        on a directory BCL files to generate Fastq files.
        :param run_dir: The 'run' directory that contains BCL files.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: The path to a sample-sheet.
        :param queue_name: The name of the Torque queue to use for processing.
        :param node_count: The number of nodes to request.
        :param nprocs: The maximum number of parallel processes to use.
        :param wall_time_limit: A hard time limit to bound processing.
        :param bcl_tool_path: The path to either bcl2fastq or bcl-convert.
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        '''
        super().__init__(run_dir,
                         output_path,
                         'ConvertJob',
                         [bcl_tool_path],
                         1000,
                         modules_to_load=modules_to_load)

        self.sample_sheet_path = sample_sheet_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.nprocs = nprocs
        self.wall_time_limit = wall_time_limit
        self.pmem = pmem
        self.bcl_tool = bcl_tool_path
        self.qiita_job_id = qiita_job_id
        self.job_script_path = join(self.output_path, f"{self.job_name}.sh")

        tmp = False
        for executable_name in ['bcl2fastq', 'bcl-convert']:
            if executable_name in self.bcl_tool:
                tmp = True
                break

        if not tmp:
            raise PipelineError(f'{self.bcl_tool} is not the path to a known'
                                'executable')

        self._file_check(self.sample_sheet_path)

        klss = KLSampleSheet(self.sample_sheet_path)

        if not validate_and_scrub_sample_sheet(klss):
            # if sample sheet is valid, we can assume it has the parameters
            # that we will need. No need to check for individual params.
            raise PipelineError(
                f"Sample sheet {self.sample_sheet_path} is not valid.")

        self._generate_job_script()

    def _generate_job_script(self):
        '''
        Generate a Torque job script for processing supplied root_directory.
        :return: The path to the newly-created job-script.
        '''
        lines = []

        lines.append("#!/bin/bash")
        lines.append(f"#PBS -N {self.qiita_job_id}_{self.job_name}")
        lines.append(f"#PBS -q {self.queue_name}")
        lines.append("#PBS -l nodes=%d:ppn=%d" % (self.node_count,
                                                  self.nprocs))
        lines.append("#PBS -V")
        lines.append("#PBS -l walltime=%d:00:00" % self.wall_time_limit)

        # send an email to the list of users defined below when a job starts,
        # terminates, or aborts. This is used to confirm that the package's
        # own reporting mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # list of users to be contacted independently of this package's
        # notification system, when a job starts, terminates, or gets aborted.
        lines.append("#PBS -M qiita.help@gmail.com")

        lines.append(f"#PBS -l pmem={self.pmem}")

        lines.append("set -x")
        lines.append('date')
        lines.append('hostname')
        lines.append(f'cd {self.root_dir}')

        if self.modules_to_load:
            lines.append("module load " + ' '.join(self.modules_to_load))

        # Assume that the bcl-convert tool is named 'bcl-convert' and choose
        # accordingly.
        if 'bcl-convert' in self.bcl_tool:
            lines.append(('%s '
                          '--sample-sheet %s '
                          '--output-directory %s '
                          '--bcl-input-directory . '
                          '--bcl-num-decompression-threads 16 '
                          '--bcl-num-conversion-threads 16 '
                          '--bcl-num-compression-threads 16 '
                          '--bcl-num-parallel-tiles 16 '
                          '--bcl-sampleproject-subdirectories true '
                          '--force') % (self.bcl_tool,
                                        self.sample_sheet_path,
                                        self.output_path))

            # equivalent cp for bcl-conversion (see below) needed.
        else:
            lines.append(('%s '
                          '--sample-sheet %s '
                          '--minimum-trimmed-read-length 1 '
                          '--mask-short-adapter-reads 1 '
                          '-R . '
                          '-o %s '
                          '--loading-threads 16 '
                          '--processing-threads 16 '
                          '--writing-threads 16 '
                          '--create-fastq-for-index-reads '
                          '--ignore-missing-positions ') %
                         (self.bcl_tool,
                          self.sample_sheet_path,
                          self.output_path))

        with open(self.job_script_path, 'w') as f:
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub(r'\s+', ' ', line)
                f.write(f"{line}\n")

    def run(self, callback=None):
        '''
        Run BCL2Fastq/BCLConvert conversion
        :param callback: optional function taking two parameters (id, status)
                         that is called when a running process's status is
                         changed.
        :return:
        '''
        job_info = self.qsub(self.job_script_path, None, None,
                             exec_from=self.log_path, callback=callback)

        logging.info(f'Successful job: {job_info}')

    def audit(self, sample_ids):
        '''
        Audit the results of a run.
        :param sample_ids: A list of sample-ids that require results.
        :return: A list of sample-ids that were not found.
        '''
        files_found = []

        for root, dirs, files in walk(self.output_path):
            files_found += [join(root, x) for x in files if
                            x.endswith('fastq.gz')]

        # remove all files found with a 'zero_files' directory in the path
        files_found = [x for x in files_found if 'zero_files' not in x]

        found = []
        for sample_id in sample_ids:
            for found_file in files_found:
                # the trailing underscore is important as it can be assumed
                # that all fastq.gz files will begin with sample_id followed
                # by an '_', and then one or more additional parameters
                # separated by underscores. This substring is unlikely to be
                if basename(found_file).startswith('%s_' % sample_id):
                    found.append(sample_id)
                    break

        return sorted(list(set(found) ^ set(sample_ids)))
