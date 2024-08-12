from jinja2 import BaseLoader, TemplateNotFound
from metapool import load_sample_sheet
from os import stat, makedirs, rename
from os.path import join, basename, dirname, exists, abspath, getmtime
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
from sequence_processing_pipeline.Pipeline import Pipeline
from shutil import move
import logging
from sequence_processing_pipeline.Commands import split_similar_size_bins
from sequence_processing_pipeline.util import iter_paired_files
from jinja2 import Environment
import glob
import re
from sys import executable
import pathlib


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




from os.path import join, exists
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
import logging
import re













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
        # being passed to TRConvertJob, additional validation isn't needed.

        self._generate_job_script()

    def _generate_script_one(self):
        template = self.jinja_env.get_template("tellread.sh")

        tellread_map = "/home/qiita_test/qiita-spots/tellread_mapping.csv"
        seqrun_path = "/sequencing/igm_runs/240216_LH00444_0058_A22357VLT4"
        lane = 'L008'
        reference_map = ""
        reference_base = ""
        mode = "metagenomic"

        return template.render(tellread_map=tellread_map,
                               seqrun_path=seqrun_path,
                               lane=lane,
                               reference_map=reference_map,
                               reference_base=reference_base,
                               mode=mode)

    def _generate_script_two(self):
        template = self.jinja_env.get_template("tellread-cleanup.sbatch")

        job_name = "cleanup"
        wall_time_limit = "24:00:00"
        mem_in_gb = "8"
        node_count = "1"
        cores_per_task = "1"
        queue_name = "qiita"

        return template.render(job_name=job_name,
                               wall_time_limit=wall_time_limit,
                               mem_in_gb=mem_in_gb,
                               node_count=node_count,
                               cores_per_task=cores_per_task,
                               queue_name=queue_name)

    def _generate_script_three(self):
        template = self.jinja_env.get_template("tellread.sbatch")
        job_name = "tellread"
        wall_time_limit = "96:00:00"
        mem_in_gb = "16"
        node_count = "1"
        cores_per_task = "4"
        queue_name = "qiita"
        tellread_sbatch_tmp_dir = "/panfs/${USER}/tmp"
        tr_sing_script_path = "$HOME/qiita-spots/tellread-release-novaseqX/run_tellread_sing.sh"
        modules_to_load = ["singularity_3.6.4"]

        return template.render(job_name=job_name,
                               wall_time_limit=wall_time_limit,
                               mem_in_gb=mem_in_gb,
                               node_count=node_count,
                               cores_per_task=cores_per_task,
                               queue_name=queue_name,
                               tmp_dir=tellread_sbatch_tmp_dir,
                               sing_script_path=tr_sing_script_path,
                               modules_to_load=' '.join(modules_to_load))

    def _generate_script_four(self):
        template = self.jinja_env.get_template("telllink-isolate.sbatch")

        job_name = "tellink-isolate"
        wall_time_limit = "96:00:00"
        node_count = "1"
        cores_per_task = "16"
        mem_in_gb = "160"
        queue_name = "qiita"
        modules_to_load = ["singularity_3.6.4"]
        sing_path = "/projects/long_read_collab/code/tellseq/release_v1.11/tellink-release/run_tellink_sing.sh"

        return template.render(job_name=job_name,
                               wall_time_limit=wall_time_limit,
                               mem_in_gb=mem_in_gb,
                               node_count=node_count,
                               cores_per_task=cores_per_task,
                               queue_name=queue_name,
                               modules_to_load=' '.join(modules_to_load),
                               sing_path=sing_path)

    def _generate_script_five(self):
        template = self.jinja_env.get_template("telllink.sbatch")

        job_name = "tellink"
        mem_in_gb = "160"
        node_count = "1"
        cores_per_task = "16"
        wall_time_limit = "96:00:00"
        queue_name = "qiita"
        modules_to_load = ["singularity_3.6.4"]
        sing_path = "/projects/long_read_collab/code/tellseq/release_v1.11/tellink-release/run_tellink_sing.sh"

        return template.render(job_name=job_name,
                               mem_in_gb=mem_in_gb,
                               node_count=node_count,
                               cores_per_task=cores_per_task,
                               wall_time_limit=wall_time_limit,
                               queue_name=queue_name,
                               modules_to_load=' '.join(modules_to_load),
                               sing_path=sing_path)

    def _generate_script_six(self):
        template = self.jinja_env.get_template("integrate.sbatch")

        job_name = "integrate"
        mem_in_gb = "8"
        node_count = "1"
        cores_per_task = "1"
        wall_time_limit = "24:00:00"
        queue_name = "qiita"

        return template.render(job_name=job_name,
                               mem_in_gb=mem_in_gb,
                               node_count=node_count,
                               cores_per_task=cores_per_task,
                               wall_time_limit=wall_time_limit,
                               queue_name=queue_name)

    def _generate_script_seven(self):
        template = self.jinja_env.get_template("cloudspades-isolate.sbatch")

        job_name = "cs-assemble"
        mem_in_gb = "64"
        node_count = "1"
        cores_per_task = "12"
        wall_time_limit = "24:00:00"
        queue_name = "qiita"
        modules_to_load = ["gcc_9.3.0"]
        spades_path = "~/spades-cloudspades-paper/assembler/spades.py"

        return template.render(job_name=job_name,
                               mem_in_gb=mem_in_gb,
                               node_count=node_count,
                               cores_per_task=cores_per_task,
                               wall_time_limit=wall_time_limit,
                               queue_name=queue_name,
                               modules_to_load=' '.join(modules_to_load),
                               spades_path=spades_path)

    def _generate_script_eight(self):
        template = self.jinja_env.get_template("cloudspades.sbatch")

        job_name = "cs-assemble"
        wall_time_limit = "24:00:00"
        mem_in_gb = "128"
        node_count = "1"
        cores_per_task = "12"
        queue_name = "qiita"
        modules_to_load = ["gcc_9.3.0"]
        spades_path = "TBD" # for now pass but don't use this spades_path var.

        return template.render(job_name=job_name,
                               mem_in_gb=mem_in_gb,
                               node_count=node_count,
                               cores_per_task=cores_per_task,
                               wall_time_limit=wall_time_limit,
                               queue_name=queue_name,
                               modules_to_load=' '.join(modules_to_load),
                               spades_path=spades_path)

    def _generate_job_scripts(self):
        scripts = [
            {
                "template": self.jinja_env.get_template("cloudspades.sbatch"),
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
            {},
            {}

        ]

        for script in scripts:
            template = self.jinja_env.get_template(script["template"])
            params = script["params"]
            result = template.render(**params)






    def _generate_job_script(self):
        """
        Generate a Torque job script for processing supplied root_directory.
        :return: The path to the newly-created job-script.
        """
        lines = []

        lines.append("#!/bin/bash")
        lines.append(f"#SBATCH --job-name {self.qiita_job_id}_{self.job_name}")
        lines.append(f"#SBATCH -p {self.queue_name}")
        lines.append(f'#SBATCH -N {self.node_count}')
        lines.append(f'#SBATCH -n {self.nprocs}')
        lines.append("#SBATCH --time %d" % self.wall_time_limit)

        # send an email to the list of users defined below when a job starts,
        # terminates, or aborts. This is used to confirm that the package's
        # own reporting mechanism is reporting correctly.
        lines.append("#SBATCH --mail-type=ALL")

        # list of users to be contacted independently of this package's
        # notification system, when a job starts, terminates, or gets aborted.
        lines.append("#SBATCH --mail-user qiita.help@gmail.com")

        lines.append(f"#SBATCH --mem-per-cpu {self.pmem}")

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
                          '--sample-sheet "%s" '
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
                          '--sample-sheet "%s" '
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
