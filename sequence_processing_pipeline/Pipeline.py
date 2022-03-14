from json import load as json_load
from json.decoder import JSONDecodeError
from os import makedirs, listdir
from os.path import join, exists, isdir
from metapool import KLSampleSheet, quiet_validate_and_scrub_sample_sheet
from metapool.plate import ErrorMessage, WarningMessage
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
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

    def __init__(self, configuration_file_path, run_id, sample_sheet_path,
                 output_path, qiita_job_id, config_dict=None):
        """
        Initialize Pipeline object w/configuration information.
        :param configuration_file_path: Path to configuration.json file.
        :param run_id: Used w/search_paths to locate input run_directory.
        :param sample_sheet_path: Path to required sample-sheet.
        :param output_path: Path where all pipeline-generated files live.
        :param qiita_job_id: Qiita Job ID creating this Pipeline.
        :param config_dict: (Optional) Dict used instead of config file.
        """
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
        for key in ['search_paths', 'archive_path']:
            if key not in config:
                raise PipelineError(f"'{key}' is not a key in "
                                    f"{self.configuration_file_path}")

        # The sample-sheet parameter is now supplied to the Pipeline object
        # at initialization time. This allows the Pipeline to process the
        # sample-sheet once and retain the metadata for subsequent method
        # calls.
        if sample_sheet_path is None:
            raise PipelineError('sample_sheet_path cannot be None')

        # Note that any warnings and/or Error messages from our extended
        # validate() method will now be returned w/

        # If our extended validate() method discovers any warnings or
        # Errors, it will raise a PipelineError and return them w/in the
        # error message as a single string separated by '\n'.
        self.warnings = []
        self.sample_sheet = self._validate_sample_sheet(sample_sheet_path)
        self._directory_check(output_path, create=False)
        self.output_path = output_path
        self.search_paths = config['search_paths']
        self.run_id = run_id
        self.run_dir = self._search_for_run_dir()
        self.qiita_job_id = qiita_job_id
        self.pipeline = []

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

    def _directory_check(self, directory_path, create=False):
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

    def _validate_sample_sheet(self, sample_sheet_path):
        """
        Performs additional validation for sample-sheet on top of metapool.
        :return: If successful, an empty list of strings and a valid
                 sample-sheet. If unsuccessful, a list of warning and error
                 messages and None.
        """
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
                # return a valid sample-sheet, and preserve any warning
                # messages
                self.warnings += [str(x) for x in msgs if
                                  isinstance(x, WarningMessage)]
                return val_sheet

        # if we are here, then val_sheet is None and there are msgs or the
        # sample-sheet failed to pass our additional tests. In either case we
        # should raise a PipelineError and return the list of ErrorMessages in
        # the PipelineError's message member.
        #
        # convert msgs from a list of ErrorMessages into a list of strings
        # before raising the PipelineError. (WarningMessages are included).
        raise PipelineError('Sample-sheet has the following errors:\n'
                            '\n'.join([str(x) for x in msgs]))

    def generate_sample_information_files(self):
        """
        Generate sample-information files in self.output_path.
        :return: A list of paths to sample-information-files.
        """
        samples = []
        for sample in self.sample_sheet.samples:
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

    def get_sample_ids(self):
        return [x.Sample_ID for x in self.sample_sheet.samples
                if 'BLANK' not in x]

    def get_project_info(self):
        bioinformatics = self.sample_sheet.Bioinformatics
        results = []

        for result in bioinformatics.to_dict('records'):
            results.append({'project_name': result['Sample_Project'],
                            'qiita_id': result['QiitaID']})

        return results
