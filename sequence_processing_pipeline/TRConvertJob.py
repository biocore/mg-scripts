from jinja2 import BaseLoader, TemplateNotFound, Environment
from os.path import split, join, exists, getmtime
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
import pathlib
from os import rename, walk, chmod, listdir, makedirs
from shutil import move, rmtree
from re import match


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
        self.pmem = pmem
        self.bcl_tool = bcl_tool_path
        self.qiita_job_id = qiita_job_id
        self.suffix = 'fastq.gz'

        self.tellread_output_path = join(self.output_path, 'output')
        makedirs(self.tellread_output_path)

        self.tmp1_path = join(self.tellread_output_path, 'tmp1')

        makedirs(self.tmp1_path)

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

        # TODO: generate a sample-mapping to map C#s to fake sample-names and
        #  fake projects. Process sample-sheet later.
        self.mapping = self._generate_sample_mapping()

        # TODO: hardcode lane at 'L001'
        self.lane = 'L001'

        self.clean_wall_time_limit = "24:00:00"
        self.clean_mem_in_gb = "8"
        self.clean_node_count = "1"
        self.clean_cores_per_task = "1"
        self.cloudspades_cores_per_task = "12"
        self.cloudspades_mem_in_gb = "128"
        self.cloudspades_modules = ["gcc_9.3.0"]
        self.cloudspades_node_count = "1"
        self.cloudspades_path = ("/home/qiita_test/qiita-spots/spades-"
                                 "cloudspades-0.1")
        self.cloudspades_wall_time_limit = "24:00:00"
        self.counts_cores_per_task = "1"
        self.counts_create_picklist_path = ("/home/qiita_test/qiita-spots/"
                                            "create_picklist.py")
        self.counts_mem_in_gb = "8"
        self.counts_node_count = "1"
        self.counts_other_file = ('20230906_FS10001773_68_BTR67708-1611.'
                                  'read_counts.tsv')
        self.counts_plot_counts_path = ("/home/qiita_test/qiita-spots/'"
                                        "'plot_counts.py")
        self.counts_sample_sheet = ("/home/qiita_test/qiita-spots/"
                                    "20230906_FS10001773_68_BTR67708-1611.csv")
        self.counts_wall_time_limit = "24:00:00"
        self.cs_isolate_mem_in_gb = "64"
        self.integrate_indicies_script_path = ("/home/qiita_test/qiita-spots/"
                                               "integrate-indices-np.py")
        self.integrate_mem_in_gb = "8"
        self.integrate_node_count = "1"
        self.integrate_wall_time_limit = "24:00:00"
        self.integrate_cores_per_task = "1"
        self.queue_name = "qiita"
        self.tellink_cores_per_task = "16"
        self.tellink_mem_in_gb = "160"
        self.tellink_modules = ["singularity_3.6.4"]
        self.tellink_node_count = "1"
        self.tellink_sing_path = ("/projects/long_read_collab/code/tellseq/"
                                  "release_v1.11/tellink-release/"
                                  "run_tellink_sing.sh")
        self.tellink_wall_time_limit = "96:00:00"
        self.tellread_cores_per_task = "4"
        self.tellread_mem_in_gb = "16"
        self.tellread_modules = ["singularity_3.6.4"]
        self.tellread_node_count = "1"
        self.tellread_sing_script_path = ("$HOME/qiita-spots/tellread-release"
                                          "-novaseqX/run_tellread_sing.sh")
        self.tellread_wall_time_limit = "96:00:00"
        self.tl_cores_per_task = "16"
        self.tl_isolate_node_count = "1"
        self.tl_isolate_wall_time_limit = "96:00:00"
        self.tl_mem_in_gb = "160"
        self.main_map = ("/home/qiita_test/qiita-spots/20230906_FS10001773_"
                         "68_BTR67708-1611.csv")
        self.main_mode = "metagenomic"
        self.main_seqrun_path = ("/sequencing/seqmount/KL_iSeq_Runs/20230906"
                                 "_FS10001773_68_BTR67708-1611")

        # TODO: Address reference_map and reference_base
        self.main_reference_base = ""
        self.main_reference_map = ""

        self._generate_job_scripts()

    def _generate_job_scripts(self):
        scripts = [
            {
                "template": "cloudspades.sbatch",
                "params": {
                    "job_name": "cs-assemble",
                    "wall_time_limit": self.wall_time_limit,
                    "mem_in_gb": self.cloudspades_mem_in_gb,
                    "node_count": self.cloudspades_node_count,
                    "cores_per_task": self.cloudspades_cores_per_task,
                    "queue_name": self.queue_name,
                    "modules_to_load": ' '.join(self.cloudspades_modules),
                    "cloudspades_path": self.cloudspades_path
                }
            },
            {
                "template": "cloudspades-isolate.sbatch",
                "params": {
                    "job_name": "cs-assemble",
                    "wall_time_limit": self.cloudspades_wall_time_limit,
                    "mem_in_gb": self.cs_isolate_mem_in_gb,
                    "node_count": self.cloudspades_node_count,
                    "cores_per_task": self.cloudspades_cores_per_task,
                    "queue_name": self.queue_name,
                    "modules_to_load": ' '.join(self.cloudspades_modules),
                    "cloudspades_path": self.cloudspades_path
                }
            },
            {
                "template": "integrate.sbatch",
                "params": {
                    "job_name": "integrate",
                    "wall_time_limit": self.integrate_wall_time_limit,
                    "mem_in_gb": self.integrate_mem_in_gb,
                    "node_count": self.integrate_node_count,
                    "cores_per_task": self.integrate_cores_per_task,
                    "iinp_script_path": self.integrate_indicies_script_path,
                    "queue_name": self.queue_name
                }
            },
            {
                "template": "compute_sequence_counts_for_normalization.sbatch",
                "params": {
                    "job_name": "norm",
                    "wall_time_limit": self.counts_wall_time_limit,
                    "mem_in_gb": self.counts_mem_in_gb,
                    "node_count": self.counts_node_count,
                    "cores_per_task": self.counts_cores_per_task,
                    "sample_sheet": self.counts_sample_sheet,
                    "plot_counts_path": self.counts_plot_counts_path,
                    "output_path": self.tellread_output_path,
                    "create_picklist_path": self.counts_create_picklist_path,
                    "read_counts_path": join(self.tellread_output_path,
                                             self.counts_other_file),
                    "queue_name": self.queue_name
                }
            },
            {
                "template": "telllink.sbatch",
                "params": {
                    "job_name": "tellink",
                    "wall_time_limit": self.tellink_wall_time_limit,
                    "mem_in_gb": self.tellink_mem_in_gb,
                    "node_count": self.tellink_node_count,
                    "cores_per_task": self.tellink_cores_per_task,
                    "queue_name": self.queue_name,
                    "modules_to_load": ' '.join(self.tellink_modules),
                    "output_path": self.tellread_output_path,
                    "sing_path": self.tellink_sing_path
                }
            },
            {
                "template": "telllink-isolate.sbatch",
                "params": {
                    "job_name": "tellink-isolate",
                    "wall_time_limit": self.tellink_wall_time_limit,
                    "node_count": self.tl_isolate_node_count,
                    "cores_per_task": self.tl_cores_per_task,
                    "mem_in_gb": self.tl_mem_in_gb,
                    "queue_name": self.queue_name,
                    "modules_to_load": ' '.join(self.tellink_modules),
                    "output_path": self.tellread_output_path,
                    "sing_path": self.tellink_sing_path
                }
            },
            {
                "template": "tellread.sbatch",
                "params": {
                    "job_name": "tellread",
                    "wall_time_limit": self.tellread_wall_time_limit,
                    "mem_in_gb": self.tellread_mem_in_gb,
                    "node_count": self.tellread_node_count,
                    "tmp_dir": self.tmp1_path,
                    "cores_per_task": self.tellread_cores_per_task,
                    "queue_name": self.queue_name,
                    "sing_script_path": self.tellread_sing_script_path,
                    "modules_to_load": ' '.join(self.tellread_modules)
                }
             },
            {
                "template": "tellread-cleanup.sbatch",
                "params": {
                    "job_name": "cleanup",
                    "wall_time_limit": self.clean_wall_time_limit,
                    "mem_in_gb": self.clean_mem_in_gb,
                    "node_count": self.clean_node_count,
                    "cores_per_task": self.clean_cores_per_task,
                    "queue_name": self.queue_name
                }
             },
            # these hardcoded paths for tellread.sh need to be replaced with
            # the lane number and run-directory path, and the lane and the
            # mode from the user input. Note that we also need to process the
            # upcoming sample-sheet in order to generate the mapping we need
            # as well.
            {
                "template": "tellread.sh",
                "params": {
                    "tellread_map": self.main_map,
                    "seqrun_path": self.main_seqrun_path,
                    "output_path": self.tellread_output_path,
                    "lane": self.lane,
                    "reference_map": self.main_reference_map,
                    "reference_base": self.main_reference_base,
                    "mode": self.main_mode
                }
             }
        ]

        for script in scripts:
            template = self.jinja_env.get_template(script["template"])
            params = script["params"]
            job_script_path = join(self.output_path, script["template"])

            with open(job_script_path, 'w') as f:
                f.write(template.render(**params))
                # TODO: Change from 777 to something more appropriate.
                chmod(job_script_path, 0o777)

    def run(self, callback=None):
        """
        Run BCL2Fastq/BCLConvert conversion
        :param callback: optional function taking two parameters (id, status)
                         that is called when a running process's status is
                         changed.
        :return:
        """

        # Unlike other Jobs that submit a Slurm script and wait for the job
        # to complete, this Job() will execute an existing shell script that
        # spawns all the jobs that perform the actual work.

        # tellread.sh performs some work that requires it to run on a compute
        # node. Since Job()s run on the interactive node, an interactive
        # shell on a compute node must be requested for this script to run on.

        # define 'sjob' here for clarity. This should be more than adequate
        # resources to run the tellread.sh script and exit as it does not wait
        # on its children to complete.

        # as with the original scripts, the scripts generated by Jinja2 will
        # live in the current working directory. Hence, the script will always
        # exist at ./tellread.sh provided it was created successfully.
        sjob = "srun -N 1 -n 1 -p qiita --mem 4g --time 1:00:00 --pty bash -l"
        command = (f"{sjob}; pushd .;cd {self.output_path}; ./tellread.sh; "
                   "popd; exit")

        if not exists(join(self.output_path, 'tellread.sh')):
            raise PipelineError("tellread.sh script could not be found.")

        res = self._system_call(command)

        if res['return_code'] != 0:
            raise PipelineError("tellread.sh script did not execute correctly")

        # once _system_call() returns and tellread.sh executed correctly, then
        # a pids file should exist in the output subdirectory.
        pids_fp = join(self.output_path, 'output', 'pids')
        if not exists(pids_fp):
            raise PipelineError("TRConvertJob could not locate a pids file")

        with open(pids_fp, 'r') as f:
            lines = f.readlines()
            lines = [x.strip().split(': ') for x in lines]
            results = {k: v for (k, v) in lines}

        child_processes = [('main tellread', 'TRJOB_RETURN_CODE',
                            'TRJOB_PID', True),
                           ('counts', 'NORM_COUNTS_JOB_RETURN_CODE',
                            'NORM_COUNTS_JOB_PID', False),
                           ('integrate', 'INTEGRATE_JOB_RETURN_CODE',
                            'INTEGRATE_JOB_PID', True),
                           ('csj', 'CSJ_JOB_RETURN_CODE',
                            'CSJ_JOB_PID', False),
                           ('tlj', 'TLJ_JOB_RETURN_CODE',
                            'TLJ_JOB_PID', False),
                           ('cleanup', 'CLEANUP_JOB_RETURN_CODE',
                            'CLEANUP_JOB_PID', True)]

        # Iterate through all the TellRead script's known child processes.
        # Some children will be optional depending on the parameters given,
        # while others are required. The Job() should immediately raise an
        # error if any child (optional or not) exits unsuccessfully, however.
        for name, code, _, is_required in child_processes:
            if code in results:
                if results[code] != '0':
                    raise PipelineError(f"An error ({results[code]}) occurred "
                                        f"running {name} subprocess")
            else:
                if is_required:
                    raise PipelineError(f"The {name} subprocess did not "
                                        "execute correctly")

        # Get a list of Slurm job ids that we need to wait on and text
        # descriptions of what they are.
        jids = [(results[x[2]], x[0]) for x in child_processes if
                x[2] in results]

        # ensure the jids are casted to integers before passing them.
        statuses = self.wait_on_job_ids([int(x[0]) for x in jids])

        for jid, description in jids:
            status = statuses[jid]
            if status not in Job.slurm_status_successful:
                raise PipelineError(f"process '{description}' ({jid}) "
                                    f"failed ({status})")

        # post-process working directory to make it appear like results
        # generated by ConvertJob

        integrated_files_path = join(self.output_path, 'output', "integrated")

        if not exists(integrated_files_path):
            raise ValueError(f"{integrated_files_path} does not exist")

        # move integrated directory to TRConvertJob directory, co-level with
        # output directory. This makes it easier to delete the rest of the
        # output that we don't need.

        # move err and out logs into logs subdirectory.
        for root, dirs, files in walk(self.output_path):
            for _file in files:
                _path = join(root, _file)
                if _path.endswith('.err'):
                    move(_path, join(self.output_path, 'logs'))
                elif _path.endswith('.out'):
                    move(_path, join(self.output_path, 'logs'))
            # don't go below one level.
            break

        # save two logs and move them into standard Job logs directory.
        move(join(self.output_path, 'output', 'log'),
             join(self.output_path, 'logs'))
        move(join(self.output_path, 'output', 'output.log'),
             join(self.output_path, 'logs'))

        # rename the files and move them into project directories.
        for root, dirs, files in walk(integrated_files_path):
            for _file in files:
                fastq_file = join(root, _file)
                self._post_process_file(fastq_file, self.mapping)

        # move project folders from integrated directory to working_dir.
        contents = listdir(integrated_files_path)
        for name in contents:
            move(join(integrated_files_path, name),
                 self.output_path)

        # delete the original output directory.
        rmtree(join(self.output_path, 'output'))

    def parse_logs(self):
        raise PipelineError("parsing logs not implemented.")

    @staticmethod
    def parse_job_script(job_script_path):
        raise PipelineError("parsing job script not implemented.")

    def _post_process_file(self, fastq_file, mapping):
        # generate names of the form generated by bcl-convert/bcl2fastq:
        # <Sample_ID>_S#_L00#_<R# or I#>_001.fastq.gz
        # see:
        # https://help.basespace.illumina.com/files-used-by-basespace/
        # fastq-files
        _dir, _file = split(fastq_file)

        # ex: integrated/C544.R2.fastq.gz
        m = match(r"(C5\d\d)\.([R,I]\d)\.fastq.gz", _file)

        if m is None:
            raise ValueError(f"The filename '{_file}' is not of a "
                             "recognizable form")

        adapter_id = m[1]
        read_type = m[2]

        if adapter_id not in mapping:
            raise ValueError(f"{adapter_id} is not present in mapping")

        sample_name, sample_index, project_name = mapping[adapter_id]

        # generate the new filename for the fastq file, and reorganize the
        # files by project.
        new_name = "%s_S%d_%s_%s_001.fastq.gz" % (sample_name,
                                                  sample_index,
                                                  self.lane,
                                                  read_type)

        # ensure that the project directory exists before we rename and move
        # the file to that location.
        makedirs(join(_dir, project_name), exist_ok=True)

        # if there's an error renaming and moving the file, let it pass up to
        # the user.
        final_path = join(_dir, project_name, new_name)
        rename(fastq_file, final_path)
        return final_path

    def _generate_sample_mapping(self):
        # this generates a sample mapping for the C501-C596 adapters used by
        # the vendor to a sample-name and project. In production use this
        # mapping would need to be created from the future sample-sheet.
        project_names = ['Project1', 'Project2', 'Project3']
        sample_mapping = {}

        for sample_index in range(1, 97):
            adapter_id = "C%s" % str(sample_index + 500)
            sample_name = "MySample%d" % sample_index
            project_name = project_names[sample_index % 3]
            sample_mapping[adapter_id] = (sample_name, sample_index,
                                          project_name)

        return sample_mapping
