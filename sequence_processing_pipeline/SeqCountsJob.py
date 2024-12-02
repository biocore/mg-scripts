from os.path import join, split
from .Job import Job, KISSLoader
from .PipelineError import JobFailedError
import logging
from jinja2 import Environment
from os import walk
from json import dumps


logging.basicConfig(level=logging.DEBUG)


class SeqCountsJob(Job):
    def __init__(self, run_dir, output_path, queue_name,
                 node_count, wall_time_limit, jmem, modules_to_load,
                 qiita_job_id, max_array_length, files_to_count_path,
                 cores_per_task=4):
        """
        ConvertJob provides a convenient way to run bcl-convert or bcl2fastq
        on a directory BCL files to generate Fastq files.
        :param run_dir: The 'run' directory that contains BCL files.
        :param output_path: Path where all pipeline-generated files live.
        :param queue_name: The name of the Torque queue to use for processing.
        :param node_count: The number of nodes to request.
        :param wall_time_limit: A hard time limit (in min) to bound processing.
        :param jmem: String representing total memory limit for entire job.
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        :param max_array_length: A hard-limit for array-sizes
        :param files_to_count_path: A path to a list of file-paths to count.
        :param cores_per_task: (Optional) # of CPU cores per node to request.
        """
        super().__init__(run_dir,
                         output_path,
                         'SeqCountsJob',
                         [],
                         max_array_length,
                         modules_to_load=modules_to_load)

        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.cores_per_task = cores_per_task

        # raise an Error if jmem is not a valid floating point value.
        self.jmem = str(int(jmem))
        self.qiita_job_id = qiita_job_id
        self.jinja_env = Environment(loader=KISSLoader('templates'))

        self.job_name = (f"seq_counts_{self.qiita_job_id}")
        self.files_to_count_path = files_to_count_path

        with open(self.files_to_count_path, 'r') as f:
            lines = f.readlines()
            lines = [x.strip() for x in lines]
            lines = [x for x in lines if x != '']
            self.file_count = len(lines)

    def run(self, callback=None):
        job_script_path = self._generate_job_script()
        params = ['--parsable',
                  f'-J {self.job_name}',
                  f'--array 1-{self.sample_count}']
        try:
            self.job_info = self.submit_job(job_script_path,
                                            job_parameters=' '.join(params),
                                            exec_from=None,
                                            callback=callback)

            logging.debug(f'SeqCountsJob Job Info: {self.job_info}')
        except JobFailedError as e:
            # When a job has failed, parse the logs generated by this specific
            # job to return a more descriptive message to the user.
            info = self.parse_logs()
            # prepend just the message component of the Error.
            info.insert(0, str(e))
            raise JobFailedError('\n'.join(info))

        self._aggregate_counts()

        logging.debug(f'SeqCountJob {self.job_info["job_id"]} completed')

    def _generate_job_script(self):
        job_script_path = join(self.output_path, "seq_counts.sbatch")
        template = self.jinja_env.get_template("seq_counts.sbatch")

        #  got to make files_to_count.txt and put it in the output directory

        with open(job_script_path, mode="w", encoding="utf-8") as f:
            f.write(template.render({
                "job_name": "seq_counts",
                "wall_time_limit": self.wall_time_limit,
                "mem_in_gb": self.jmem,
                "node_count": self.node_count,
                "cores_per_task": self.cores_per_task,
                "queue_name": self.queue_name,
                "file_count": self.file_count,
                "output_path": self.output_path
            }))

        return job_script_path

    def parse_logs(self):
        # TODO
        pass

    def _aggregate_counts(self):
        def extract_metadata(fp):
            with open(fp, 'r') as f:
                lines = f.readlines()
                lines = [x.strip() for x in lines]
                if len(lines) != 2:
                    raise ValueError("error processing %s" % fp)
                _dir, _file = split(lines[0])
                seq_counts, base_pairs = lines[1].split('\t')
                return _dir, _file, int(seq_counts), int(base_pairs)

        results = {}

        for root, dirs, files in walk(self.log_path):
            for _file in files:
                if _file.endswith('.out'):
                    log_output_file = join(root, _file)
                    _dir, _file, seq_counts, base_pairs = \
                        extract_metadata(log_output_file)

                    if _dir not in results:
                        results[_dir] = {}

                    results[_dir][_file] = {'seq_counts': seq_counts,
                                            'base_pairs': base_pairs}

        results_path = join(self.output_path, 'aggregate_counts.json')

        with open(results_path, 'w') as f:
            print(dumps(results, indent=2), file=f)

        return results_path
