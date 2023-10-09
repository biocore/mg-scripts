from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os import stat, listdir, makedirs
from os.path import join
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
from shutil import move
import logging
from datetime import date
from sequence_processing_pipeline.Commands import split_similar_size_bins
from jinja2 import Environment, FileSystemLoader


logging.basicConfig(level=logging.DEBUG)


'''
TODO:
The only two methods not called in testing are:
in _filter
in run
'''


class NuQCJob(Job):
    def __init__(self, fastq_root_dir, output_path, sample_sheet_path,
                 minimap_database_paths, queue_name,
                 node_count, nprocs, wall_time_limit, jmem, fastp_path,
                 minimap2_path, samtools_path, modules_to_load, qiita_job_id,
                 pool_size, max_array_length, known_adapters_path):
        """
        Submit a Torque job where the contents of fastq_root_dir are processed
        using fastp, minimap, and samtools. Human-genome sequences will be
        filtered out, if needed.
        :param fastq_root_dir: Path to a dir of Fastq files, org. by project.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: Path to a sample sheet file.
        :param minimap_database_paths: Path to human genome databases in env.
        :param queue_name: Torque queue name to use in running env.
        :param node_count: Number of nodes to use in running env.
        :param nprocs: Number of processes to use in runing env.
        :param wall_time_limit: Hard wall-clock-time limit (in min) for procs.
        :param jmem: String representing total memory limit for entire job.
        :param fastp_path: The path to the fastp executable
        :param minimap2_path: The path to the minimap2 executable
        :param samtools_path: The path to the samtools executable
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        :param pool_size: The number of jobs to process concurrently.
        :param known_adapters_path: The path to an .fna file of known adapters.

        """
        super().__init__(fastq_root_dir,
                         output_path,
                         'QCJob',
                         [fastp_path, minimap2_path, samtools_path],
                         max_array_length,
                         modules_to_load=modules_to_load)
        self.sample_sheet_path = sample_sheet_path
        self._file_check(self.sample_sheet_path)
        metadata = self._process_sample_sheet()
        self.sample_ids = metadata['sample_ids']
        self.project_data = metadata['projects']
        self.needs_trimming = metadata['needs_adapter_trimming']
        self.nprocs = 16 if nprocs > 16 else nprocs
        self.chemistry = metadata['chemistry']
        # instead of a list of individual .mmi files,  minimap_database_paths
        # is now a path to a directory of .mmi files.
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
        self.jinja_env = Environment(loader=FileSystemLoader("templates/"))

        # 'fastp_known_adapters_formatted.fna'
        self.known_adapters_path = known_adapters_path

        ### BEGIN NEW PARAMS
        # sizes in GB
        # max size of GB data to process in a job
        self.max_file_list_size_in_gb = 2

        datestamp = date.today().strftime("%Y.%m.%d")

        self.batch_prefix = join(self.output_path,
                                 f"1-hd-split-pangenome-pe-{datestamp}")

        self.temp_dir = join(self.output_path, 'tmp')
        makedirs(self.temp_dir, exist_ok=True)
        ### END NEW PARAMS

        self.minimum_bytes = 3100

        if not isinstance(self.needs_trimming, bool):
            raise ValueError("needs_adapter_trimming must be boolean.")

        # These values are no longer needed as Daniel and Caitlin's new
        # method handles these situations more intelligently. However this
        # code is needed to identify and raise Errors when sample-sheets and
        # pre-prep files are given invalid values.

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

        # TODO: runprefix could contain _R1_ and so this should be rwritten to
        #  be a little more robust. Also there can be instances of .R1. in
        #  file-names.
        for entry in listdir(filtered_directory):
            if '_R1_' in entry:
                reverse_entry = entry.replace('_R1_', '_R2_')
                full_path = join(filtered_directory, entry)
                full_path_reverse = join(filtered_directory, reverse_entry)
                if stat(full_path).st_size <= minimum_bytes or stat(
                        full_path_reverse).st_size <= minimum_bytes:
                    logging.debug(f'moving {entry} and {reverse_entry}'
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

        batch_count = split_similar_size_bins(self.root_dir,
                                              self.max_file_list_size_in_gb,
                                              self.batch_prefix)

        export_params = [f"MMI={self.minimap_database_paths}",
                         f"PREFIX={self.batch_prefix}",
                         f"OUTPUT={self.output_path}",
                         f"TMPDIR={self.temp_dir}"]

        job_params = ['-J', self.batch_prefix, '--array', batch_count,
                      '--mem', '20G', '--export', ','.join(export_params)]

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
        logging.debug(f'QCJob {job_id} completed')

        '''
        for project in self.project_data:
            project_name = project['Sample_Project']
            needs_human_filtering = project['HumanFiltering']

            source_dir = join(self.output_path, project_name)


            # TODO: UNCONFIRMED COMPLETED FILES WILL BE CREATED
            # confirm Slurm job was successful by using .completed files.
            # if not self._get_failed_indexes(project_name, job_id):
            #     raise PipelineError("QCJob did not complete successfully.")

            ### DANIEL'S NEW CODE DOES NOT APPEAR TO FILTER FOR
                ZERO-LENGTH FILES ###

            # determine where the filtered fastq files can be found and move
            # the 'zero-length' files to another directory.
            if needs_human_filtering is True:
                filtered_directory = join(source_dir, 'filtered_sequences')
            else:
                filtered_directory = join(source_dir, 'trimmed_sequences')

            empty_files_directory = join(source_dir, 'zero_files')
            self._filter_empty_fastq_files(filtered_directory,
                                           empty_files_directory,
                                           self.minimum_bytes)
        '''

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

        # Adapter trimming is not appropriate for metatranscriptomics, so it
        # is not included here.
        needs_adapter_trimming = header['Assay'] == Pipeline.METAGENOMIC_PTYPE

        sample_ids = []
        for sample in valid_sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        bioinformatics = valid_sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.
        lst = bioinformatics.to_dict('records')

        # convert true/false and yes/no strings to true boolean values.
        for record in lst:
            for key in record:
                if record[key].strip().lower() in ['true', 'yes']:
                    record[key] = True
                elif record[key].strip().lower() in ['false', 'no']:
                    record[key] = False

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry,
                'projects': lst,
                'sample_ids': sample_ids,
                'needs_adapter_trimming': needs_adapter_trimming
                }

    def _generate_job_script(self):
        job_script_path = join(self.output_path, 'process_all_fastq_files.sh')
        template = self.jinja_env.get_template("nuqc_job.sh")

        job_name = f'{self.qiita_job_id}_{self.job_name}'

        with open(job_script_path, mode="w", encoding="utf-8") as f:
            f.write(template.render(job_name=job_name,
                                    wall_time_limit_in_min=5760,
                                    # for NuQC, mem should be 20GB
                                    mem_in_gb=self.jmem,
                                    # should be 1
                                    node_count=self.node_count,
                                    cores_per_task=8,
                                    known_adapters_path=self.known_adapters_path))

        return job_script_path
