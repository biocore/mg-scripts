from jinja2 import BaseLoader, TemplateNotFound, Environment
from os.path import join, exists, getmtime
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
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
        self.job_script_path = None

        # for projects that use sequence_processing_pipeline as a dependency,
        # jinja_env must be set to sequence_processing_pipeline's root path,
        # rather than the project's root path.
        self.jinja_env = Environment(loader=KISSLoader('templates'),
                                     # set Jinja2 comment strings to be
                                     # anything other than '{#' and '#}',
                                     # which can be used in shell scripts.
                                     comment_start_string='%%%%%%%%%%',
                                     comment_end_string='%%%%%%%%%%')

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
                    "job_name": "tellink",
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
                    "tmp_dir": "/panfs/${USER}/tmp",
                    "cores_per_task": "4",
                    "queue_name": "qiita",
                    "tellread_sbatch_tmp_dir": "/panfs/${USER}/tmp",
                    "sing_script_path": "$HOME/qiita-spots/tellread-release"
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
                "template": "tellread.sh",
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

            if script['template'] == "tellread.sh":
                self.job_script_path = job_script_path

    def run(self, callback=None):
        """
        Run BCL2Fastq/BCLConvert conversion
        :param callback: optional function taking two parameters (id, status)
                         that is called when a running process's status is
                         changed.
        :return:
        """

        # Unlike other Jobs that submit a Slurm script and wait for the job
        # to complete, this Job will execute the tellread.sh shell script.
        # It is this script that does all of the Slurm job creation. This Job
        # will need another means to tell when a job has completed
        # successfully.

        command = ("./tellread.sh -s /sequencing/igm_runs/240216_LH00444"
                   "_0058_A22357VLT4 -i ./samplesheet.csv -l L008 -m "
                   "metagenomic")

        if self.job_script_path:
            res = self._system_call(command)
        else:
            raise PipelineError("tellread.sh script could not be found.")

        if res['return_code'] != 0:
            raise PipelineError("tellread.sh script did not execute correctly")

        # res['stdout']
        # res['stderr']

    def parse_logs(self):
        raise PipelineError("parsing logs not implemented.")

    @staticmethod
    def parse_job_script(job_script_path):
        raise PipelineError("parsing job script not implemented.")
