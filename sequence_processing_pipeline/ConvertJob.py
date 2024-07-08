from os.path import join, exists, split
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError)
import logging
import re
from json import loads as json_loads
from metapool import load_sample_sheet
from shutil import copyfile


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
        self.pmem = pmem
        self.bcl_tool = bcl_tool_path
        self.qiita_job_id = qiita_job_id
        self.job_script_path = join(self.output_path, f"{self.job_name}.sh")
        self.suffix = 'fastq.gz'
        self.fastq_paths = None
        self.info = None

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

            self._get_sample_sheet_info()

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

    def _get_sample_sheet_info(self):
        # assume path to sample-sheet exists, and sheet is valid.
        # otherwise, we would not be here.
        sheet = load_sample_sheet(self.sample_sheet_path)

        # parse bioinformatics section to generate a durable list of
        # project_names and qiita_ids.
        bioinformatics = sheet.Bioinformatics
        projects = bioinformatics['Sample_Project'].tolist()
        qiita_ids = bioinformatics['QiitaID'].tolist()

        if 'contains_replicates' in sheet.Bioinformatics:
            has_reps = sheet.Bioinformatics['contains_replicates'].tolist()
            # assume a validated sample-sheet ensures has_reps has only one
            # value, either True or False.
            self.contains_replicates = bool(has_reps[0])
        else:
            self.contains_replicates = False

        results = {}

        for project, qiita_id in zip(projects, qiita_ids):
            # derive project_name by removing the known qiita_id associated
            # with this project from the string.
            project_name = re.sub(f'_{qiita_id}$', '', project)
            results[project] = {'qiita_id': qiita_id,
                                'project_name': project_name,
                                'full_project_name': project,
                                'samples': {}}

        # Since the project-name is stored in an internal variable
        # in a third-party library, convert samples to JSON using the exposed
        # method first.
        samples = json_loads(sheet.to_json())['Data']

        for sample in samples:
            d = {'Sample_Name': sample['Sample_Name'],
                 'Sample_ID': sample['Sample_ID'],
                 # matching files will store the paths to all fastq files
                 # associated w/this sample-name.
                 'matching_files': []}

            if 'orig_name' in sample:
                d['orig_name'] = sample['orig_name']

            results[sample['Sample_Project']]['samples'][d['Sample_Name']] = d

        # associate with each dictionary a list of matching fastq files.
        # this way, we can easily determine which files to copy to another
        # project based on just a sample-name/sample-id and a project-name.

        for project in results:
            # find just the fastq files for this project.

            fastq_paths = self._find_files(join(self.output_path,
                                                project))
            fastq_paths = [f for f in fastq_paths if f.endswith('.fastq.gz')]

            for sample_name in results[project]['samples']:
                sample = results[project]['samples'][sample_name]
                # regex based on studying all filenames of all fastq files in
                # $WKDIR. Works with _R1_, _R2_, _I1_, _I2_, etc.
                rgx = r"^" + re.escape(sample['Sample_ID']) + \
                      r"_S\d+_L\d+_[R,I]\d+_\d+.fastq.gz$"

                for full_path in fastq_paths:
                    file_path, file_name = split(full_path)
                    if re.match(rgx, file_name):
                        sample['matching_files'].append(full_path)

        self.info = results

    def copy_sequences(self, sample_name, source_project, dest_project,
                       copy_all_replicates=False):
        """
        Copies all fastq files related to a sample into another project.
        :param source_project: The source project w/qiita_id.
        :param dest_project: The destination project w/qiita_id.
        :param sample_name: A sample-name.
        :param orig_name: A sample-name.
        :param copy_all_replicates: If True, search for sample_name in the
            orig_name column of the sample-sheet. Copy all replicates.
        :return: None
        """
        if self.info is None:
            raise ValueError("This method cannot be called until processing "
                             "has completed.")

        project_names = list(self.info.keys())

        # confirm source project is a valid one.
        if source_project not in project_names:
            raise ValueError(f"'{source_project}' is not defined in the "
                             "sample-sheet")

        # confirm destination project is a valid one.
        if dest_project not in project_names:
            raise ValueError(f"'{dest_project}' is not defined in the "
                             "sample-sheet")

        if source_project == dest_project:
            raise ValueError(f"source '{source_project}' and destination "
                             f"'{dest_project}' projects are the same")

        # note that the user can supply a sample-name that didn't make it
        # through the conversion step and may have no files matched to it.
        # this is considered okay and a possible outcome of conversion. In
        # this case zero files are copied and no information is passed back
        # to the user.

        # projects that contain replicates must also be considered. if the
        # value for sample_name contains a well-id, then only the files
        # associated with that particular replicate should be copied. If the
        # value instead references a sample_name in the 'orig_name' column,
        # then all files associated with each replicate need to be moved.

        # in this situation, the sample_name needs be compared against all
        # orig_names in the project and the individual sample_names (w/well-
        # ids) must be discovered. Then those individual sample_names can
        # be processed.

        if copy_all_replicates is True and self.contains_replicates is False:
            raise ValueError("'treat_as_orig_name' is set to 'True' but this "
                             "sample-sheet doesn't contain replicates")

        samples = self.info[source_project]['samples']
        sample_list = []

        if copy_all_replicates:
            for key in samples:
                sample = samples[key]
                # assume orig_name is present if treat_as_orig_name is True.
                if sample_name == sample['orig_name']:
                    sample_list.append(sample)
        else:
            # sample_name is a value from the sample_name column. it may or
            # may not have a well-id appended and this sample-sheet may or
            # may not contain replicates, but in either case a single sample
            # either exists or it doesn't.
            if sample_name in self.info[source_project]['samples']:
                sample_list.append(samples[sample_name])

        if len(sample_list) == 0:
            # if the sample_list is empty, then sample-name wasn't present in
            # either the sample_name or orig_name columns.
            raise ValueError(f"'{sample_name}' is not defined in the project"
                             f" '{source_project}'")

        for sample in sample_list:
            for src_fp in sample['matching_files']:
                # split(fp)[1] is simply the original filename, which must
                # be provided in the destination path.
                dst_fp = join(self.output_path, dest_project, split(src_fp)[1])
                copyfile(src_fp, dst_fp)
