from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os import stat, makedirs
from os.path import join, basename, dirname, exists
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
from shutil import move
import logging
from sequence_processing_pipeline.Commands import split_similar_size_bins
from sequence_processing_pipeline.util import iter_paired_files
from jinja2 import Environment, PackageLoader
import glob
import re
from json import dumps
from sys import executable


logging.basicConfig(level=logging.DEBUG)


class NuQCJob(Job):
    def __init__(self, fastq_root_dir, output_path, sample_sheet_path,
                 minimap_database_paths, queue_name, node_count,
                 wall_time_limit, jmem, fastp_path, minimap2_path,
                 samtools_path, modules_to_load, qiita_job_id, pool_size,
                 max_array_length, known_adapters_path, bucket_size=8):
        """
        Submit a slurm job where the contents of fastq_root_dir are processed
        using fastp, minimap2, and samtools. Human-genome sequences will be
        filtered out, if needed.
        :param fastq_root_dir: Path to a dir of Fastq files, org. by project.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: Path to a sample sheet file.
        :param minimap_database_paths: Path to human genome databases in env.
        :param queue_name: Torque queue name to use in running env.
        :param node_count: Number of nodes to use in running env.
        :param wall_time_limit: Hard wall-clock-time limit (in min) for procs.
        :param jmem: String representing total memory limit for entire job.
        :param fastp_path: The path to the fastp executable
        :param minimap2_path: The path to the minimap2 executable
        :param samtools_path: The path to the samtools executable
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        :param pool_size: The number of jobs to process concurrently.
        :param known_adapters_path: The path to an .fna file of known adapters.
        :param bucket_size: the size in GB of each bucket to process

        """
        super().__init__(fastq_root_dir,
                         output_path,
                         'NuQCJob',
                         [fastp_path, minimap2_path, samtools_path],
                         max_array_length,
                         modules_to_load=modules_to_load)
        self.sample_sheet_path = sample_sheet_path
        self._file_check(self.sample_sheet_path)
        metadata = self._process_sample_sheet()
        self.sample_ids = metadata['sample_ids']
        self.project_data = metadata['projects']
        self.chemistry = metadata['chemistry']
        self.minimap_database_paths = minimap_database_paths
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.jmem = jmem
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size
        self.suffix = 'fastq.gz'

        # for projects that use sequence_processing_pipeline as a dependency,
        # jinja_env must be set to sequence_processing_pipeline's root path,
        # rather than the project's root path.
        self.jinja_env = Environment(loader=PackageLoader('sequence_processing'
                                                          '_pipeline',
                                                          'templates'))

        self.counts = {}
        self.known_adapters_path = known_adapters_path
        self.max_file_list_size_in_gb = bucket_size
        self.temp_dir = join(self.output_path, 'tmp')
        makedirs(self.temp_dir, exist_ok=True)
        self.batch_prefix = "hd-split-pangenome"
        self.minimum_bytes = 3100
        self._validate_project_data()

    def _validate_project_data(self):
        # Validate project settings in [Bioinformatics]
        for project in self.project_data:
            if project['ForwardAdapter'] == 'NA':
                project['ForwardAdapter'] = None

            if project['ReverseAdapter'] == 'NA':
                project['ReverseAdapter'] = None

            if project['ForwardAdapter'] is None:
                if project['ReverseAdapter'] is not None:
                    raise ValueError(("ForwardAdapter is declared but not "
                                      "ReverseAdapter."))

            if project['ReverseAdapter'] is None:
                if project['ForwardAdapter'] is not None:
                    raise ValueError(("ReverseAdapter is declared but not "
                                      "ForwardAdapter."))

            if not isinstance(project['HumanFiltering'], bool):
                raise ValueError("needs_adapter_trimming must be boolean.")

    def _filter_empty_fastq_files(self, filtered_directory,
                                  empty_files_directory,
                                  minimum_bytes):
        '''
        Filters out and moves fastq files that are below threshold.
        :param filtered_directory:
        :param empty_files_directory:
        :param minimum_bytes:
        :return:
        '''
        empty_list = []

        files = glob.glob(join(filtered_directory, f'*.{self.suffix}'))

        for r1, r2 in iter_paired_files(files):
            full_path = join(filtered_directory, r1)
            full_path_reverse = join(filtered_directory, r2)
            if stat(full_path).st_size <= minimum_bytes or stat(
                    full_path_reverse).st_size <= minimum_bytes:
                logging.debug(f'moving {full_path} and {full_path_reverse}'
                              f' to empty list.')
                empty_list.append(full_path)
                empty_list.append(full_path_reverse)

        if empty_list:
            logging.debug(f'making directory {empty_files_directory}')
            makedirs(empty_files_directory, exist_ok=True)

        for item in empty_list:
            logging.debug(f'moving {item}')
            move(item, empty_files_directory)

    def run(self, callback=None):
        # now a single job-script will be created to process all projects at
        # the same time, and intelligently handle adapter-trimming as needed
        # as well as human-filtering.
        job_script_path = self._generate_job_script()

        batch_location = join(self.temp_dir, self.batch_prefix)
        batch_count = split_similar_size_bins(self.root_dir,
                                              self.max_file_list_size_in_gb,
                                              batch_location)
        self.counts[self.batch_prefix] = batch_count

        export_params = [f"MMI={self.minimap_database_paths}",
                         f"PREFIX={self.batch_prefix}",
                         f"OUTPUT={self.output_path}",
                         f"TMPDIR={self.temp_dir}"]

        job_params = ['-J', self.batch_prefix, f'--array 1-{batch_count}',
                      '--export', ','.join(export_params)]

        # job_script_path formerly known as:
        #  process.multiprep.pangenome.adapter-filter.pe.sbatch
        job_info = self.submit_job(job_script_path,
                                   # job_parameters - was None
                                   ' '.join(job_params),
                                   # script_parameters
                                   None,
                                   # assume we want to exec from the log_path
                                   # for now.
                                   exec_from=self.log_path,
                                   callback=callback)
        job_id = job_info['job_id']
        logging.debug(f'NuQCJob {job_id} completed')

        for project in self.project_data:
            project_name = project['Sample_Project']
            needs_human_filtering = project['HumanFiltering']

            source_dir = join(self.output_path, project_name)
            pattern = f"{source_dir}/*.fastq.gz"
            completed_files = list(glob.glob(pattern))

            if needs_human_filtering is True:
                filtered_directory = join(source_dir, 'filtered_sequences')
            else:
                filtered_directory = join(source_dir, 'trimmed_sequences')

            # create the properly named directory to move files to in
            # in order to preserve legacy behavior.
            makedirs(filtered_directory, exist_ok=True)

            # get the list of sample-names in this project.
            samples_in_project = [x[0] for x in self.sample_ids
                                  if x[1] == project_name]

            files_to_move = []
            # Tissue_1_Mag_Hom_DNASe_RIBO_S16_L001_R2_001.fastq.gz
            # Nislux_SLC_Trizol_DNASe_S7_L001_R2_001.fastq.gz
            regex = re.compile(r'^(.*)_S\d{1,2}_L\d{3}_R\d_\d{3}\.fastq\.gz$')
            for fp in completed_files:
                file_name = basename(fp)
                substr = regex.search(file_name)
                if substr is None:
                    raise ValueError(f"{file_name} does not follow naming "
                                     " pattern.")
                else:
                    # check if found substring is a member of this
                    # project. Note sample-name != sample-id
                    if substr[1] in samples_in_project:
                        files_to_move.append(fp)

            for fp in files_to_move:
                move(fp, filtered_directory)

            # once fastq.gz files have been moved into the right project,
            # we now need to consider the html and json fastp_reports
            # files.
            old_html_path = join(self.output_path, 'fastp_reports_dir', 'html')
            old_json_path = join(self.output_path, 'fastp_reports_dir', 'json')

            new_html_path = join(source_dir, 'fastp_reports_dir', 'html')
            new_json_path = join(source_dir, 'fastp_reports_dir', 'json')

            makedirs(new_html_path, exist_ok=True)
            makedirs(new_json_path, exist_ok=True)

            pattern = f"{old_html_path}/*.html"
            completed_htmls = list(glob.glob(pattern))

            for fp in completed_htmls:
                move(fp, new_html_path)

            pattern = f"{old_json_path}/*.json"
            completed_jsons = list(glob.glob(pattern))

            for fp in completed_jsons:
                move(fp, new_json_path)

            # now that files are separated by project as per legacy
            # operation, continue normal processing.
            empty_files_directory = join(source_dir, 'zero_files')
            self._filter_empty_fastq_files(filtered_directory,
                                           empty_files_directory,
                                           self.minimum_bytes)

    def _get_failed_indexes(self, project_name, job_id):
        pattern = f"{self.temp_dir}/{self.batch_prefix}.*.completed"
        completed_files = list(glob.glob(pattern))
        completed_indexes = []
        regex = r'^%s/%s.%s_([0-9]+).completed$' % (self.temp_dir,
                                                    self.batch_prefix,
                                                    str(job_id))
        array_ids = re.compile(regex)

        for completed_file in completed_files:
            match = array_ids.search(completed_file)
            if match is None:
                raise PipelineError("Malformed completed file")
            else:
                id_ = int(match.groups(0)[0])
                completed_indexes.append(id_)

        # a successfully completed job array should have a list of array
        # numbers from 0 - len(self.commands).
        all_indexes = list(range(1, self.counts[self.batch_prefix]))
        failed_indexes = sorted(set(all_indexes) - set(completed_indexes))

        # generate log-file here instead of in run() where it can be
        # unittested more easily.
        log_fp = join(self.output_path,
                      'logs',
                      f'failed_indexes_{job_id}.json')

        if failed_indexes:
            with open(log_fp, 'w') as f:
                f.write(dumps({'job_id': job_id,
                               'failed_indexes': failed_indexes}, indent=2))

        return failed_indexes

    def _process_sample_sheet(self):
        sheet = KLSampleSheet(self.sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % self.sample_sheet_path
            raise PipelineError(s)

        header = valid_sheet.Header
        chemistry = header['chemistry']

        if header['Assay'] not in Pipeline.assay_types:
            s = "Assay value '%s' is not recognized." % header['Assay']
            raise PipelineError(s)

        sample_ids = []
        for sample in valid_sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        bioinformatics = valid_sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.
        lst = bioinformatics.to_dict('records')

        # convert true/false and yes/no strings to true boolean values.
        for record in lst:
            # the subset of columns that should be either True or False.
            for key in ['BarcodesAreRC', 'HumanFiltering']:
                val = record[key].strip()
                if val == 'True':
                    record[key] = True
                elif val == 'False':
                    record[key] = False
                else:
                    raise ValueError(f"'{val}' is not a valid value for {key}")

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry,
                'projects': lst,
                'sample_ids': sample_ids}

    def _generate_job_script(self):
        job_script_path = join(self.output_path, 'process_all_fastq_files.sh')
        template = self.jinja_env.get_template("nuqc_job.sh")

        job_name = f'{self.qiita_job_id}_{self.job_name}'

        html_path = join(self.output_path, 'fastp_reports_dir', 'html')
        json_path = join(self.output_path, 'fastp_reports_dir', 'json')

        # get location of python executable in this environment.
        # demux script should be present in the same location.
        demux_path = join(dirname(executable), 'demux')

        if not exists(demux_path):
            raise ValueError(f"{demux_path} does not exist.")

        with open(job_script_path, mode="w", encoding="utf-8") as f:
            # the job resources should come from a configuration file
            f.write(template.render(job_name=job_name,
                                    queue_name=self.queue_name,
                                    # should be 4 * 24 * 60 = 4 days
                                    wall_time_limit=self.wall_time_limit,
                                    mem_in_gb=self.jmem,
                                    node_count=1,
                                    cores_per_task=4,
                                    knwn_adpt_path=self.known_adapters_path,
                                    output_path=self.output_path,
                                    html_path=html_path,
                                    json_path=json_path,
                                    demux_path=demux_path))

        return job_script_path
