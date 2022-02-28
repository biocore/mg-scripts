from json import load as json_load
from json.decoder import JSONDecodeError
from os import makedirs, listdir, walk
from os.path import basename, join, exists, isdir, getmtime
from metapool import KLSampleSheet, quiet_validate_and_scrub_sample_sheet
from metapool.plate import ErrorMessage
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from time import time as epoch_time
import logging
import re


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)


class Pipeline:
    sif_header = ['sample_name', 'collection_timestamp', 'elevation', 'empo_1',
                  'empo_2', 'empo_3', 'env_biome', 'env_feature',
                  'env_material', 'env_package', 'geo_loc_name',
                  'host_subject_id', 'latitude', 'longitude', 'sample_type',
                  'scientific_name', 'taxon_id', 'description', 'title',
                  'dna_extracted', 'physical_specimen_location',
                  'physical_specimen_remaining']

    sif_defaults = [None, None, 193, 'Control', 'Negative',
                    'Sterile water blank', 'urban biome', 'research facility',
                    'sterile water', 'misc environment', 'USA:CA:San Diego',
                    None, 32.5, -117.25, 'control blank', 'metagenome', 256318,
                    None, 'adaptation', 'TRUE', 'UCSD', 'FALSE']

    def __init__(self, configuration_file_path, run_id, output_path,
                 qiita_job_id, config_dict=None):
        '''
        Pipeline
        :param configuration_file_path: Path to configuration.json file.
        :param run_id: Used w/search_paths to locate input run_directory.
        :param output_path: Path where all pipeline-generated files live.
        :param qiita_job_id: Qiita Job ID creating this Pipeline.
        :param config_dict: (Optional) Dict used instead of config file.
        '''
        if config_dict:
            if 'configuration' in config_dict:
                self.configuration = config_dict['configuration']
                self.configuration_file_path = None
            else:
                raise PipelineError(f"{config_dict} does not contain the "
                                    "key 'configuration'")
        else:
            self.configuration_file_path = configuration_file_path
            try:
                f = open(configuration_file_path)
                self.configuration = json_load(f)['configuration']
                f.close()
            except TypeError:
                raise PipelineError('configuration_file_path cannot be None')
            except FileNotFoundError:
                raise PipelineError(f'{configuration_file_path} does not '
                                    'exist.')
            except JSONDecodeError:
                raise PipelineError(f'{configuration_file_path} is not a '
                                    'valid json file')

        if run_id is None:
            raise PipelineError('run_id cannot be None')

        if 'pipeline' not in self.configuration:
            raise PipelineError("'pipeline' is not a key in "
                                f"{self.configuration_file_path}")

        config = self.configuration['pipeline']
        for key in ['search_paths', 'archive_path', 'younger_than',
                    'older_than']:
            if key not in config:
                raise PipelineError(f"'{key}' is not a key in "
                                    f"{self.configuration_file_path}")

        self.output_path = output_path
        self.directory_check(output_path, create=False)
        self.search_paths = config['search_paths']
        self.run_id = run_id
        self.run_dir = self._search_for_run_dir()
        self.younger_than = config['younger_than']
        self.older_than = config['older_than']
        self.run_id = run_id
        self.qiita_job_id = qiita_job_id
        self.pipeline = []

        if self.older_than < 0 or self.younger_than < 0:
            raise PipelineError('older_than and younger_than cannot be '
                                'less than zero.')

        if self.older_than >= self.younger_than:
            raise PipelineError(
                'older_than cannot be equal to or less than younger_than.')

        # required files for successful operation
        # both RTAComplete.txt and RunInfo.xml should reside in the root of
        # the run directory.
        required_files = ['RTAComplete.txt', 'RunInfo.xml']
        for some_file in required_files:
            if not exists(join(self.run_dir, some_file)):
                raise PipelineError("required file '%s' is not present." %
                                    some_file)

        # verify that RunInfo.xml file is readable.
        try:
            fp = open(join(self.run_dir, 'RunInfo.xml'))
            fp.close()
        except PermissionError:
            raise PipelineError('RunInfo.xml is present, but not readable')

    def _search_for_run_dir(self):
        # this method will catch a run directory as well as its products
        # directory, which also has the same name. Hence, return the
        # shortest matching path as that will at least return the right
        # path between the two.
        results = []

        for search_path in self.search_paths:
            logging.debug(f'Searching {search_path} for {self.run_id}')
            for entry in listdir(search_path):
                some_path = join(search_path, entry)
                # ensure some_path never ends in '/'
                some_path = some_path.rstrip('/')
                if isdir(some_path) and some_path.endswith(self.run_id):
                    logging.debug(f'Found {some_path}')
                    results.append(some_path)

        if results:
            results.sort(key=lambda s: len(s))
            return results[0]

        raise PipelineError(f"A run-dir for '{self.run_id}' could not be "
                            "found")

    def directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PipelineError, so it gets handled in the same location
                    # as the others.
                    raise PipelineError(str(e))
            else:
                raise PipelineError("directory_path '%s' does not exist." %
                                    directory_path)

    def is_within_time_range(self, run_directory):
        """
        Verify run_directory's timestamp is within Pipeline's window.
        :return: True if run_directory is w/in window, False otherwise.
        """
        if exists(run_directory):
            ts = getmtime(run_directory)
            d_t = epoch_time() - ts
            if d_t < (self.younger_than * 3600):
                if d_t > (self.older_than * 3600):
                    return True
        else:
            raise PipelineError(f"'{run_directory}' doesn't exist")

        return False

    def run(self):
        """
        Run all jobs added to Pipeline in the order they were added.
        :return: None
        """
        for job in self.pipeline:
            job.run()

    def add(self, job):
        """
        Add a job to the Pipeline
        :param Job: A Job object
        :return: None
        """
        if isinstance(job, Job):
            self.pipeline.append(job)
        else:
            raise PipelineError("object is not a Job object.")

    def validate(self, sample_sheet_path):
        '''
        Performs additional validation for sample-sheet on top of metapool.
        :param sample_sheet_path: Path to sample-sheet.
        :return: If successful, an empty list of strings and a valid
                 sample-sheet. If unsuccessful, a list of warning and error
                 messages and None.
        '''
        # validate the sample-sheet using metapool package.
        sheet = KLSampleSheet(sample_sheet_path)
        msgs, val_sheet = quiet_validate_and_scrub_sample_sheet(sheet)

        passes_additional_tests = True

        if val_sheet is not None:
            # perform extended validation based on required fields for
            # seqpro, and other issues encountered.
            bioinformatics = val_sheet.Bioinformatics
            if 'library_construction_protocol' not in bioinformatics:
                msgs.append(ErrorMessage("column 'library_construction_protoco"
                                         "l' not found in Bioinformatics secti"
                                         "on"))
            if 'experiment_design_description' not in bioinformatics:
                msgs.append(ErrorMessage("column 'experiment_design_descriptio"
                                         "n' not found in Bioinformatics secti"
                                         "on"))

            # look for duplicate samples. metapool will allow two rows w/the
            # same lane and sample_id if one or more other columns are
            # different. However seqpro expects the tuple (lane, sample_id) to
            # be unique for indexing.
            unique_indexes = []
            for item in val_sheet.samples:
                unique_index = f'{item.lane}_{item.sample_id}'
                if unique_index in unique_indexes:
                    passes_additional_tests = False
                    msgs.append(ErrorMessage("A sample already exists with la"
                                             f"ne {item.lane} and sample-id "
                                             f"{item.sample_id}"))
                else:
                    unique_indexes.append(unique_index)

        if passes_additional_tests:
            return msgs, val_sheet
        else:
            return msgs, None

    def generate_sample_information_files(self, sample_sheet_path):
        '''
        Using a path to a validated sample-sheet, generate sample-information
        files in self.output_path.
        :param sample_sheet_path:
        :return: A list of paths to sample-information-files.
        '''
        sheet = KLSampleSheet(sample_sheet_path)
        msgs, val_sheet = quiet_validate_and_scrub_sample_sheet(sheet)

        if val_sheet is None:
            raise PipelineError("'%s' is not a valid sample-sheet." %
                                sample_sheet_path)

        samples = []
        for sample in val_sheet.samples:
            if sample['Sample_ID'].startswith('BLANK'):
                samples.append((sample['Sample_ID'], sample['Sample_Project']))

        projects = list(set([y for x, y in samples]))

        paths = []
        for project in projects:
            samples_in_proj = [x for x, y in samples if y == project]
            some_path = join(self.output_path, f'{project}_blanks.tsv')
            paths.append(some_path)
            with open(some_path, 'w') as f:
                # write out header to disk
                f.write('\t'.join(Pipeline.sif_header) + '\n')

                # for now, populate values that can't be derived from the
                # sample-sheet w/'EMPTY'.
                for sample in samples_in_proj:
                    row = {}
                    for column, default_value in zip(Pipeline.sif_header,
                                                     Pipeline.sif_defaults):
                        # ensure all defaults are converted to strings.
                        row[column] = str(default_value)

                    # overwrite default title w/sample_project name, minus
                    # Qiita ID.
                    row['title'] = re.sub(r'_\d+$', r'', project)

                    # generate values for the four columns that must be
                    # determined from sample-sheet information.

                    # convert 'BLANK14_10F' to 'BLANK14.10F', etc.
                    row['sample_name'] = sample.replace('_', '.')
                    row['host_subject_id'] = sample.replace('_', '.')
                    row['description'] = sample.replace('_', '.')

                    # generate collection_timestamp from self.run_id
                    # assume all run_ids begin with coded datestamp:
                    # 210518_...
                    # allow exception if substrings cannot convert to int
                    # or if array indexes are out of bounds.
                    year = int(self.run_id[0:2]) + 2000
                    month = int(self.run_id[2:4])
                    day = int(self.run_id[4:6])
                    row['collection_timestamp'] = f'{year}-{month}-{day}'

                    row = [row[x] for x in Pipeline.sif_header]
                    f.write('\t'.join(row) + '\n')

        return paths

    def extract_metadata(self, sample_sheet_path):
        '''
        Extract sample ids from a given sample-sheet.
        :param sample_sheet_path: Path to a sample-sheet.
        :return: A dictionary of various lists of metadata.
        '''
        sheet = KLSampleSheet(sample_sheet_path)
        msgs, val_sheet = quiet_validate_and_scrub_sample_sheet(sheet)

        if val_sheet is None:
            raise PipelineError("'%s' is not a valid sample-sheet." %
                                sample_sheet_path)

        # generate a list of unique project, plate names, and sample-ids found
        # in the sample-sheet.
        projects = list(set([x.Sample_Project for x in val_sheet.samples]))
        plates = list(set([x.Sample_Plate for x in val_sheet.samples]))
        sample_ids = list(set([x.Sample_ID for x in val_sheet.samples]))

        # Generate a list of BLANKs defined on each plate. They will be shared
        # across all projects on that plate. For now, accept '' as a valid
        # plate-name.
        blanks = [x for x in sample_ids if 'BLANK' in x]

        # filter out BLANKs from the list of sample-ids.
        sample_ids = [x for x in sample_ids if 'BLANK' not in x]

        # map blanks to plates
        bp_map = {plate: [] for plate in plates}

        # map samples to projects
        spro_map = {project: [] for project in projects}

        # map samples to plates
        splate_map = {plate: [] for plate in plates}

        for sample in val_sheet.samples:
            if sample.Sample_ID in blanks:
                bp_map[sample.Sample_Plate].append(sample.Sample_ID)
            else:
                spro_map[sample.Sample_Project].append(sample.Sample_ID)
                splate_map[sample.Sample_Plate].append(sample.Sample_ID)

        result = {
            'unique-values': {
                'project-names': projects,
                'plate-names': plates,
                'sample-ids': sample_ids
            }, 'maps': {
                'blanks-to-plates': bp_map,
                'samples-to-plates': splate_map,
                'samples-to-projects': spro_map
            }
        }

        return result

    def find_files_ending_in(self, working_path, suffix):
        '''
        Find files ending in 'suffix' in a nested working directory.
        :param working_path: Root path of nested directory.
        :param suffix: Sub-string to filter for.
        :return:
        '''
        results = []

        for root, dirs, files in walk(working_path):
            results += [join(root, x) for x in files if x.endswith(suffix)]

        return results

    def audit_job(self, sample_ids, working_path, suffix, mod):
        '''
        Audit the results of a Job object.
        :param sample_ids: A list of sample-ids that require results.
        :param working_path: Root path of the Job object's working directory.
        :param suffix: The file-type or extension of the files to filter for.
        :param mod: The substring marking the end of a sample-id.
        :return: A list of sample-ids found.
        '''
        files_found = self.find_files_ending_in(working_path, suffix)
        # remove all files found with a 'zero_files' directory in the path
        files_found = [x for x in files_found if 'zero_files' not in x]
        found = []

        for sample_id in sample_ids:
            for found_file in files_found:
                # we can append '_R\d' to sample_id if we
                # can assume _R1 or _R2 will always follow the
                # sample_id
                mod_id = '%s_%s' % (sample_id, mod)
                if basename(found_file).startswith(mod_id):
                    found.append(sample_id)
                    break

        return found


if __name__ == '__main__':
    logging.debug("Starting Pipeline")

    sample_sheet_path = '/path/to/sample-sheet'
    run_id = '210518_A00953_0305_test11'
    config_file_path = '/path/to/configuration.json'
    qiita_job_id = 'NOT_A_QIITA_JOB_ID'
    output_path = '/path/to/output_path'
    config_dict = None

    pipeline = Pipeline(config_file_path, run_id, output_path, qiita_job_id,
                        config_dict)

    msgs, val_sheet = pipeline.validate(sample_sheet_path)
    paths = pipeline.generate_sample_information_files(sample_sheet_path)

    sdo = SequenceDirectory(pipeline.run_dir,
                            sample_sheet_path=sample_sheet_path)

    config = pipeline.configuration['bcl-convert']

    pipeline.add(ConvertJob(pipeline.run_dir,
                            pipeline.output_path,
                            sdo.sample_sheet_path,
                            config['queue'],
                            config['nodes'],
                            config['nprocs'],
                            config['wallclock_time_in_hours'],
                            config['per_process_memory_limit'],
                            config['executable_path'],
                            config['modules_to_load'],
                            pipeline.qiita_job_id))

    config = pipeline.configuration['qc']

    '''
    the path to the raw fastq_files will be in run_dir/ConvertJob
    '''

    raw_fastq_files_path = join(pipeline.output_path, 'ConvertJob')

    pipeline.add(QCJob(raw_fastq_files_path,
                       pipeline.output_path,
                       sdo.sample_sheet_path,
                       config['mmi_db'],
                       config['queue'],
                       config['nodes'],
                       config['nprocs'],
                       config['wallclock_time_in_hours'],
                       config['job_total_memory_limit'],
                       config['fastp_executable_path'],
                       config['minimap2_executable_path'],
                       config['samtools_executable_path'],
                       config['modules_to_load'],
                       qiita_job_id,
                       config['job_pool_size'],
                       config['job_max_array_length']))

    config = pipeline.configuration['fastqc']

    processed_fastq_files_path = join(pipeline.output_path, 'QCJob')

    # TODO: run_dir as it is isn't used by FastQCJob.
    #  Inputs are raw_fastq_files_path and processed_fastq_files_path.
    pipeline.add(FastQCJob(pipeline.run_dir,
                           pipeline.output_path,
                           raw_fastq_files_path,
                           processed_fastq_files_path,
                           config['nprocs'],
                           config['nthreads'],
                           config['fastqc_executable_path'],
                           config['modules_to_load'],
                           qiita_job_id,
                           config['queue'],
                           config['nodes'],
                           config['wallclock_time_in_hours'],
                           config['job_total_memory_limit'],
                           config['job_pool_size'],
                           config['multiqc_config_file_path'],
                           config['job_max_array_length']))

    config = pipeline.configuration['seqpro']

    pipeline.add(GenPrepFileJob(pipeline.run_dir,
                                pipeline.output_path,
                                sdo.sample_sheet_path,
                                config['seqpro_path'],
                                config['modules_to_load'],
                                qiita_job_id))

    pipeline.run()
