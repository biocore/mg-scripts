from sequence_processing_pipeline.Job import Job
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join, dirname, split
from os import walk, remove, stat, listdir, makedirs, rmdir
import logging
from sequence_processing_pipeline.QCHelper import QCHelper
from shutil import move


class QCJob(Job):
    def __init__(self, run_dir, sample_sheet_path, mmi_db_path,
                 output_directory, queue_name, node_count, nprocs,
                 wall_time_limit, pmem, fastp_path, minimap2_path,
                 samtools_path, modules_to_load):
        '''
        Submit a Torque job where the contents of run_dir are processed using
        fastp, minimap, and samtools. Human-genome sequences will be filtered
        out if needed.
        :param run_dir: Path to a run directory.
        :param sample_sheet_path: Path to a sample sheet file.
        :param fpmmp_path: Path to fpmmp.sh script in running environment.
        :param mmi_db_path: Path to human genome database in running env.
        :param output_directory: Completed runs root path. E.g.:
         /sequencing/ucsd_2/complete_runs
        :param queue_name: Torque queue name to use in running env.
        :param node_count: Number of nodes to use in running env.
        :param nprocs: Number of processes to use in runing env.
        :param wall_time_limit: Hard wall-clock-time limit for processes.
        :param pmem: String representing memory limit per process E.g.: '10gb'
        :param fastp_path: The path to the fastp executable
        :param minimap2_path: The path to the minimap2 executable
        :param samtools_path: The path to the samtools executable
        :param modules_to_load: A list of Linux module names to load
        '''
        # for now, keep this run_dir instead of abspath(run_dir)
        self.job_name = 'QCJob'
        super().__init__(run_dir,
                         self.job_name,
                         [fastp_path, minimap2_path, samtools_path],
                         modules_to_load)
        metadata = self._process_sample_sheet(sample_sheet_path)
        self.sample_ids = metadata['sample_ids']
        self.project_data = metadata['projects']
        self.needs_a_trimming = metadata['needs_adapter_trimming']
        self.trim_file = 'split_file_'
        self.nprocs = nprocs
        self.chemistry = metadata['chemistry']
        self.mmi_db_path = mmi_db_path
        self.queue_name = queue_name
        self.node_count = node_count
        self.nprocs = nprocs
        self.wall_time_limit = wall_time_limit
        self.products_dir = join(self.run_dir, 'products')
        self.pmem = pmem
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        self.modules_to_load = modules_to_load

        # POST-PROCESSING RELATED
        # self.destination_directory = '/pscratch/seq_test/test_copy'
        self.destination_directory = output_directory

        # set to 500 bytes to avoid empty and small files that Qiita
        # has trouble with.
        self.minimum_bytes = 500

    def _copy(self, reports_directory, source_directory,
              destination_directory):
        # move the fastp html and json directories up-front.
        logging.debug('moving html directory')
        move(join(reports_directory, 'html'), source_directory)
        logging.debug('moving json directory')
        move(join(reports_directory, 'json'), source_directory)

        # delete the reports_directory
        # this is because the directory should be empty, but often it is
        # also a subdirectory of the directory we will be later copying
        # to the final location.
        try:
            logging.debug('removing fastp_report_dir directory')
            rmdir(reports_directory)
        except OSError as e:
            # This will usually be because there are still directories
            # in the directory.
            raise PipelineError(str(e))

        # Now that the disk layout of the results is in the form
        # expected by users in /sequencing... Move the entire results over
        # to /sequencing or another final location w/one rsync call.
        # For performance and reliability reasons, use rsync for copying.
        self._system_call('rsync -avp %s %s' % (
            source_directory, destination_directory))

    def _filter(self, filtered_directory, empty_files_directory,
                minimum_bytes):
        empty_list = []

        for entry in listdir(filtered_directory):
            logging.debug(f'{entry} found.')
            if '_R1_' in entry:
                logging.debug(f'checking size of {entry}.')
                reverse_entry = entry.replace('_R1_', '_R2_')
                full_path = join(filtered_directory, entry)
                full_path_reverse = join(filtered_directory, reverse_entry)
                size1 = stat(full_path).st_size
                logging.debug(f'checking size of {reverse_entry}.')
                size2 = stat(full_path_reverse).st_size
                if size1 <= minimum_bytes or size2 <= minimum_bytes:
                    logging.debug(f'moving {entry} and {reverse_entry}'
                                  f' to empty list.')
                    empty_list.append(full_path)
                    empty_list.append(full_path_reverse)

        if empty_list:
            # for now, we won't mind if the empty_files_directory exists
            # already. We'll assume nothing is in it or if there are files
            # they can be overwritten and don't hurt anything.
            logging.debug(f'making directory {empty_files_directory}')
            makedirs(empty_files_directory, exist_ok=True)

        for item in empty_list:
            logging.debug(f'moving {item}')
            move(item, empty_files_directory)

    def run(self):
        for project in self.project_data:
            fastq_files = self._find_fastq_files(project['Sample_Project'])
            split_count = self._generate_split_count(len(fastq_files))
            self._clear_trim_files()
            lines_per_split = int((
                len(fastq_files) + split_count - 1) / split_count)
            trim_files = self._generate_trim_files(
                fastq_files, lines_per_split)
            if len(trim_files) != split_count:
                logging.warning("The number of trim files does not equal the "
                                "number of files we were supposed to have.")

            script_path = self._generate_job_script(self.queue_name,
                                                    self.node_count,
                                                    self.nprocs,
                                                    self.wall_time_limit,
                                                    split_count,
                                                    project['Sample_Project'],
                                                    project['ForwardAdapter'],
                                                    project['ReverseAdapter'],
                                                    self.needs_a_trimming,
                                                    project['HumanFiltering'],
                                                    self.pmem)

            pbs_job_id = self.qsub(script_path, None, None)
            logging.debug(f'QCJob {pbs_job_id} completed')

            source_dir = join(self.products_dir, project['Sample_Project'])

            # convention is for this directory to be under source_directory
            filtered_directory = join(source_dir, 'filtered_sequences')
            reports_directory = join(source_dir, 'fastp_reports_dir')

            # typically named 'zero_files'
            empty_files_directory = join(source_dir, 'zero_files')

            # perform needed filtering and copying of data once
            # fastp/minimap2/samtools have completed.
            logging.debug('possible copy to '
                          '/qmounts/qiita_data/uploads/${qiita_proj}/ '
                          "under qiita's ownership disabled.")
            self._filter(filtered_directory, empty_files_directory,
                         self.minimum_bytes)
            self._copy(reports_directory, source_dir,
                       self.destination_directory)

    def _generate_trim_files(self, fastq_files, split_count):
        def _chunk_list(some_list, n):
            # taken from https://bit.ly/38S9O8Z
            for i in range(0, len(some_list), n):
                yield some_list[i:i + n]

        # put these split files in the same location as the fastq_files for
        # the project. Assume all filepaths in fastq_files have the same
        # result for dirname().
        destination_path = dirname(fastq_files[0])

        new_files = []
        count = 0
        for chunk in _chunk_list(fastq_files, split_count):
            trim_file_name = '%s%d' % (self.trim_file, count)
            trim_file_path = join(destination_path, trim_file_name)
            with open(trim_file_path, 'w') as f:
                for line in chunk:
                    f.write("%s\n" % line)
            new_files.append(trim_file_path)
            count += 1

        return new_files

    def _clear_trim_files(self):
        # remove all files with a name beginning in self.trim_file.
        # assume cleaning the entire run_dir is overkill, but won't
        # hurt anything.
        for root, dirs, files in walk(self.run_dir):
            for some_file in files:
                if self.trim_file in some_file:
                    some_path = join(root, some_file)
                    remove(some_path)

    def _process_sample_sheet(self, sample_sheet_path):
        sheet = KLSampleSheet(sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % sample_sheet_path
            raise PipelineError(s)

        header = valid_sheet.Header
        chemistry = header['chemistry']
        needs_adapter_trimming = ('TRUE' if
                                  header['Assay'] == 'Metagenomics'
                                  else 'FALSE')

        sample_ids = []
        for sample in valid_sheet.samples:
            logging.debug(f"Found sample: {sample['Sample_ID']} in project"
                          f"{sample['Sample_Project']}.")
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        bioinformatics = valid_sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.
        lst = bioinformatics.to_dict('records')

        # convert true/false and yes/no strings to true boolean values.
        for record in lst:
            for key in ['BarcodesAreRC', 'HumanFiltering']:
                if record[key].strip().lower() in ['true', 'yes']:
                    record[key] = True
                elif record[key].strip().lower() in ['false', 'no']:
                    record[key] = False
                else:
                    raise PipelineError(f"Unexpected value '{record[key]}' "
                                        f"for column '{key}' in sample-sheet '"
                                        f"{sample_sheet_path}'")

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry,
                'projects': lst,
                'sample_ids': sample_ids,
                'needs_adapter_trimming': needs_adapter_trimming
                }

    def _find_fastq_files_in_run_dir(self, project_name):
        search_path = join(self.run_dir, 'Data', 'Fastq', project_name)
        lst = []
        for root, dirs, files in walk(search_path):
            for some_file in files:
                some_file = some_file.decode('utf-8')
                if some_file.endswith('fastq.gz'):
                    some_path = join(search_path, some_file)
                    lst.append(some_path)
        return lst

    def _find_fastq_files(self, project_name):
        # filter the list of (sample_id, sample_project) tuples stored in
        # self.sample_ids so that only the ids matching project_name are in
        # the list.
        sample_ids = filter(lambda c: c[1] == project_name, self.sample_ids)
        # strip out the project name from the matching elements.
        sample_ids = [x[0] for x in sample_ids]

        # The sample-sheet only contains sample IDs, not actual filenames.
        # Hence, we need to generate a list of possible fastq files to process,
        # and filter out the ones that don't contain samples mentioned in the
        # file.
        files_found = self._find_fastq_files_in_run_dir(project_name)
        logging.debug(f"Fastq files found: {files_found}")

        lst = []

        for some_file in files_found:
            # use os.path.split() to reliably return the file's name.
            file_path, file_name = split(some_file)
            # assume each file_name begins with a sample_id. Note that this
            # name is generated by programs like bcl-convert and bcl2fastq
            # so this is fairly durable.
            # e.g.: X00185910 and X00185910_S177_L003_R1_001.fastq.gz.
            tmp = file_name.split('_')[0]
            if tmp in sample_ids:
                # this Fastq file is one we should process.
                logging.debug(f"Located Fastq file for id '{id}': {some_file}")
                lst.append(some_file)   # Save the full path.
        return lst

    def _generate_split_count(self, count):
        if count > 2000:
            return 16
        elif count <= 2000 and count > 1000:
            return 10
        elif count <= 1000 and count > 500:
            return 4

        return 1

    def _generate_job_script(self, queue_name, node_count, nprocs,
                             wall_time_limit, split_count, project_name,
                             adapter_a, adapter_A, a_trim, h_filter, pmem):
        lines = []

        # unlike w/ConvertBCL2FastqJob, multiple qc.sh scripts will be
        # generated, one for each project defined in the sample sheet.
        # since we will generate a job-script for each project, we won't use
        # Job.stdout_log_path and stderr_log_path. Instead we'll use the
        # numbered ones provided by generate_job_script_path().
        tmp1, tmp2, tmp3 = self.generate_job_script_path()
        job_script_path = tmp1
        output_log_path = tmp2
        error_log_path = tmp3

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying
        # params w/environment variables is to do the same on QSUB but with
        # -v instead of --export. The syntax of the values are otherwise the
        # same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        # lines.append("#PBS -N {}".format(project_name))
        lines.append("#PBS -N %s" % project_name)

        # what torque calls a queue, slurm calls a partition
        # SBATCH -p SOMETHING -> PBS -q SOMETHING
        # (Torque doesn't appear to have a quality of service (-q) option. so
        # this will go unused in translation.)
        lines.append("#PBS -q %s" % queue_name)

        # request one node
        # Slurm --ntasks-per-node=<count> -> -l ppn=<count>	in Torque
        lines.append("#PBS -l nodes=%d:ppn=%d" % (node_count, nprocs))

        # Slurm --export=ALL -> Torque's -V
        lines.append("#PBS -V")

        # Slurm
        # walltime limit --time=24:00:00 -> Torque's -l walltime=<hh:mm:ss>
        # using the larger value found in the two scripts (72 vs ? hours)
        lines.append("#PBS -l walltime=%d:00:00" % wall_time_limit)

        # send an email to the list of users defined below when a job starts,
        # terminates, or aborts. This is used to confirm that the package's
        # own reporting mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # list of users to be contacted independently of this package's
        # notification system, when a job starts, terminates, or gets aborted.
        lines.append("#PBS -M "
                     "ccowart@ucsd.edu,"
                     "jdereus@ucsd.edu,"
                     "qiita.help@gmail.com")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        lines.append(f"#PBS -l pmem={pmem}")

        # --output -> -o
        lines.append("#PBS -o %s" % output_log_path)
        lines.append("#PBS -e %s" % error_log_path)

        # array input files are labeled file0, file1,...filen-1
        lines.append("#PBS -t 0-%d" % (split_count - 1))

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1

        # We probably do not need to activate this Python environment, but I
        # will store it here in comments.
        # source ~/miniconda3/bin/activate test_env_2

        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use run_dir instead.
        lines.append("set -x")
        lines.append("cd %s" % self.run_dir)
        tmp = "module load " + ' '.join(self.modules_to_load)
        logging.debug(f"QCJob Modules to load: {tmp}")
        lines.append(tmp)

        project_products_dir = join(self.products_dir, project_name)

        qc = QCHelper(self.nprocs, self.trim_file, project_name,
                      project_products_dir, self.mmi_db_path, adapter_a,
                      adapter_A, a_trim, h_filter, self.chemistry,
                      self.fastp_path, self.minimap2_path, self.samtools_path)

        # append the commands generated by fpmmp to the job script.
        lines += qc.generate_commands()

        with open(job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                f.write("%s\n" % line)

        return job_script_path
