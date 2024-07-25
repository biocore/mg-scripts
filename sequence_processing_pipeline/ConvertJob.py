from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
import logging
import re
from jinja2 import BaseLoader, Environment, TemplateNotFound
import pathlib
from os.path import join, exists, getmtime


# taken from https://jinja.palletsprojects.com/en/3.0.x/api/#jinja2.BaseLoader
class KISSLoader(BaseLoader):
    def __init__(self, path):
        # pin the path for loader to the location sequence_processing_pipeline
        # (the location of this file), along w/the relative path to the
        # templates directory.
        self.path = join(pathlib.Path(__file__).parent.resolve(), path)

    def get_source(self, environment, template):
        path = join(self.path, template)
        if not exists(path):
            raise TemplateNotFound(template)
        mtime = getmtime(path)
        with open(path) as f:
            source = f.read()
        return source, path, lambda: mtime == getmtime(path)


logging.basicConfig(level=logging.DEBUG)


class ConvertJob(Job):
    def __init__(self, run_dir, output_path, sample_sheet_path, queue_name,
                 node_count, nprocs, wall_time_limit, pmem, bcl_tool_path,
                 modules_to_load, qiita_job_id):
        """
        ConvertJob provides a convenient way to run bcl-convert or bcl2fastq
        on a directory BCL files to generate Fastq files.
        :param run_dir: The 'run' directory that contains BCL files.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: The path to a sample-sheet.
        :param queue_name: The name of the Torque queue to use for processing.
        :param node_count: The number of nodes to request.
        :param nprocs: The maximum number of parallel processes to use.
        :param wall_time_limit: A hard time limit (in min) to bound processing.
        :param bcl_tool_path: The path to either bcl2fastq or bcl-convert.
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        """
        super().__init__(run_dir,
                         output_path,
                         'ConvertJob',
                         [bcl_tool_path],
                         1000,
                         modules_to_load=modules_to_load)

        # for metagenomics pipelines, sample_sheet_path will reflect a real
        # sample_sheet file. For amplicon pipelines, sample_sheet_path will
        # reference a dummy sample_sheet file.
        self.sample_sheet_path = sample_sheet_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.nprocs = nprocs
        self.wall_time_limit = wall_time_limit

        # TODO: This value is currently a string e.g.: '1gb' or '10gb' read
        #  in from the configuration json file. However this param should be
        #  changed to process_mem_in_gb or similar and the string changed to
        #  a numerical value.
        self.pmem = pmem
        self.bcl_tool = bcl_tool_path
        self.qiita_job_id = qiita_job_id
        # CHARLIE
        self.job_script_path = join(self.output_path, f"{self.job_name}.sh")
        self.suffix = 'fastq.gz'

        # for projects that use sequence_processing_pipeline as a dependency,
        # jinja_env must be set to sequence_processing_pipeline's root path,
        # rather than the project's root path.
        self.jinja_env = Environment(loader=KISSLoader('templates'))

        tmp = False
        for executable_name in ['bcl2fastq', 'bcl-convert']:
            if executable_name in self.bcl_tool:
                tmp = True
                break

        if not tmp:
            raise PipelineError(f'{self.bcl_tool} is not the path to a known'
                                'executable')

        self._file_check(self.sample_sheet_path)

        # As the sample-sheet is validated by the Pipeline object before
        # being passed to ConvertJob, additional validation isn't needed.

        self._generate_job_script()

    def _generate_job_script(self):
        # bypass generating job script for a force-fail job, since it is
        # not needed.
        if self.force_job_fail:
            return None

        template = self.jinja_env.get_template("convert_job.sh")

        job_name = f'{self.qiita_job_id}_{self.job_name}'

        with open(self.job_script_path, mode="w", encoding="utf-8") as f:
            if 'bcl-convert' in self.bcl_tool:
                cmd_line = (f'{self.bcl_tool} '
                            f'--sample-sheet "{self.sample_sheet_path}" '
                            f'--output-directory {self.output_path} '
                            '--bcl-input-directory . '
                            '--bcl-num-decompression-threads 16 '
                            '--bcl-num-conversion-threads 16 '
                            '--bcl-num-compression-threads 16 '
                            '--bcl-num-parallel-tiles 16 '
                            '--bcl-sampleproject-subdirectories true '
                            '--force')
                # equivalent cp for bcl-conversion (see below) needed.
            else:
                cmd_line = (f'{self.bcl_tool} '
                            f'--sample-sheet "{self.sample_sheet_path}" '
                            '--minimum-trimmed-read-length 1 '
                            '--mask-short-adapter-reads 1 '
                            '-R . '
                            f'-o {self.output_path} '
                            '--loading-threads 16 '
                            '--processing-threads 16 '
                            '--writing-threads 16 '
                            '--create-fastq-for-index-reads '
                            '--ignore-missing-positions ')

            params = {'job_name': job_name,
                      'queue_name': self.queue_name,
                      'node_count': self.node_count,
                      'nprocs': self.nprocs,
                      'wall_time_limit': self.wall_time_limit,
                      'mem_per_cpu': self.pmem,
                      'run_dir': self.root_dir,
                      'sheet_path': self.sample_sheet_path,
                      'cmd_line': cmd_line}

            # generate a string of linux system modules to load before
            # processing begins.
            if self.modules_to_load:
                # if {{modules_to_load}} is defined, not empty and not false,
                # then the line "module load <modules to load>" will be
                # added to the template.
                params['modules_to_load'] = ' '.join(self.modules_to_load)

            f.write(template.render(**params))

        return self.job_script_path

    def run(self, callback=None):
        """
        Run BCL2Fastq/BCLConvert conversion
        :param callback: optional function taking two parameters (id, status)
                         that is called when a running process's status is
                         changed.
        :return:
        """
        try:
            job_info = self.submit_job(self.job_script_path,
                                       exec_from=self.log_path,
                                       callback=callback)
        except JobFailedError as e:
            # When a job has failed, parse the logs generated by this specific
            # job to return a more descriptive message to the user.
            info = self.parse_logs()
            # prepend just the message component of the Error.
            info.insert(0, str(e))
            raise JobFailedError('\n'.join(info))

        logging.info(f'Successful job: {job_info}')

    def parse_logs(self):
        log_path = join(self.output_path, 'Logs')
        errors = join(log_path, 'Errors.log')

        msgs = []

        if not exists(errors):
            # we do not raise an Error in this case because it's expected that
            # parse_logs() will be called in response to an exceptional
            # condition.
            msgs.append(f"'{errors} does not exist")

        with open(errors, 'r') as f:
            lines = f.readlines()
            for line in [x.strip() for x in lines]:
                msgs.append(line)

        return msgs

    @staticmethod
    def parse_job_script(job_script_path):
        # Returns run-directory and sample-sheet path from a job-script.

        if not exists(job_script_path):
            raise ValueError(f"'{job_script_path}' is not a valid path")

        with open(job_script_path, 'r') as f:
            lines = f.readlines()
            lines = [x.strip() for x in lines]

        # As this code creates this file, we can expect it to be of a certain
        # format.
        if lines[0] != '#!/bin/bash':
            raise ValueError(f"'{job_script_path}' is not a valid path")

        result = {}

        m = re.match('^cd (.*)$', lines[12])

        if m:
            result['run_directory'] = m.group(1)
        else:
            raise ValueError("could not detect run_directory in "
                             f"'{job_script_path}'")

        m = re.match('^bcl-convert --sample-sheet "(.*?)" ', lines[14])

        if m:
            result['sample_sheet_path'] = m.group(1)
        else:
            raise ValueError("could not detect sample-sheet path in "
                             f"'{job_script_path}'")

        return result
