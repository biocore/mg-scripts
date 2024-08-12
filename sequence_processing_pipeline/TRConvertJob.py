from jinja2 import BaseLoader, TemplateNotFound, Environment
from os.path import join, exists, getmtime
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
import logging
import pathlib
import re


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


class TRConvertJob(Job):
    def __init__(self, run_dir, output_path, sample_sheet_path, queue_name,
                 node_count, nprocs, wall_time_limit, pmem, bcl_tool_path,
                 modules_to_load, qiita_job_id):
        """
        TRConvertJob provides a convenient way to run bcl-convert or bcl2fastq
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
                         'TRConvertJob',
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
        self.pmem = pmem
        self.bcl_tool = bcl_tool_path
        self.qiita_job_id = qiita_job_id
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
        # being passed to TRConvertJob, additional validation isn't needed.

        self._generate_job_scripts()

    def _generate_job_scripts(self):
        scripts = [
            {
                "template": "cloudspades.sbatch",
                "params": {
                    "job_name": "cs-assemble",
                    "wall_time_limit": "24:00:00",
                    "mem_in_gb": "128",
                    "node_count": "1",
                    "cores_per_task": "12",
                    "queue_name": "qiita",
                    "modules_to_load": ' '.join(["gcc_9.3.0"]),
                    "spades_path": "TBD"
                }
            },
            {
                "template": "cloudspades-isolate.sbatch",
                "params": {
                    "job_name": "cs-assemble",
                    "wall_time_limit": "24:00:00",
                    "mem_in_gb": "64",
                    "node_count": "1",
                    "cores_per_task": "12",
                    "queue_name": "qiita",
                    "modules_to_load": ' '.join(["gcc_9.3.0"]),
                    "spades_path": "~/spades-cloudspades-paper/assembler/"
                                   "spades.py"
                }
            },
            {
                "template": "integrate.sbatch",
                "params": {
                    "job_name": "integrate",
                    "wall_time_limit": "24:00:00",
                    "mem_in_gb": "8",
                    "node_count": "1",
                    "cores_per_task": "1",
                    "queue_name": "qiita"
                }
            },
            {
                "template": "telllink.sbatch",
                "params": {
                    "job_name": "telllink",
                    "wall_time_limit": "96:00:00",
                    "mem_in_gb": "160",
                    "node_count": "1",
                    "cores_per_task": "16",
                    "queue_name": "qiita",
                    "modules_to_load": ' '.join(["singularity_3.6.4"]),
                    "sing_path": "/projects/long_read_collab/code/tellseq/"
                                 "release_v1.11/tellink-release/"
                                 "run_tellink_sing.sh"
                }
            },
            {
                "template": "telllink-isolate.sbatch",
                "params": {
                    "job_name": "tellink-isolate",
                    "wall_time_limit": "96:00:00",
                    "node_count": "1",
                    "cores_per_task": "16",
                    "mem_in_gb": "160",
                    "queue_name": "qiita",
                    "modules_to_load": ' '.join(["singularity_3.6.4"]),
                    "sing_path": "/projects/long_read_collab/code/tellseq/"
                                 "release_v1.11/tellink-release/"
                                 "run_tellink_sing.sh"
                }
            },
            {
                "template": "tellread.sbatch",
                "params": {
                    "job_name": "tellread",
                    "wall_time_limit": "96:00:00",
                    "mem_in_gb": "16",
                    "node_count": "1",
                    "cores_per_task": "4",
                    "queue_name": "qiita",
                    "tellread_sbatch_tmp_dir": "/panfs/${USER}/tmp",
                    "tr_sing_script_path": "$HOME/qiita-spots/tellread-release"
                                           "-novaseqX/run_tellread_sing.sh",
                    "modules_to_load": ' '.join(["singularity_3.6.4"])
                }
             },
            {
                "template": "tellread-cleanup.sbatch",
                "params": {
                    "job_name": "cleanup",
                    "wall_time_limit": "24:00:00",
                    "mem_in_gb": "8",
                    "node_count": "1",
                    "cores_per_task": "1",
                    "queue_name": "qiita"
                }
             },
            {
                "template": "",
                "params": {
                    "tellread_map": "/home/qiita_test/qiita-spots/"
                                    "tellread_mapping.csv",
                    "seqrun_path": "/sequencing/igm_runs/"
                                   "240216_LH00444_0058_A22357VLT4",
                    "lane": 'L008',
                    "reference_map": "",
                    "reference_base": "",
                    "mode": "metagenomic"
                }
             }
        ]

        for script in scripts:
            template = self.jinja_env.get_template(script["template"])
            params = script["params"]
            job_script_path = join(self.output_path, script["template"])
            with open(job_script_path, 'w') as f:
                f.write(template.render(**params))

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
