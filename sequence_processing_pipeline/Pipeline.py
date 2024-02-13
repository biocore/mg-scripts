from json import load as json_load
from json import loads as json_loads
from json.decoder import JSONDecodeError
from os import makedirs, listdir
from os.path import join, exists, isdir, basename
from metapool import load_sample_sheet, AmpliconSampleSheet
from metapool.plate import ErrorMessage, WarningMessage
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
import logging
from re import sub, findall, search
import sample_sheet
import pandas as pd
from collections import defaultdict


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)


class Pipeline:
    assay_types = ['TruSeq HT', 'Metagenomic', 'Metatranscriptomic']

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

    mapping_file_columns = {'sample_name', 'barcode', 'center_name',
                            'center_project_name',
                            'experiment_design_description',
                            'instrument_model',
                            'library_construction_protocol',
                            'platform', 'run_center', 'run_date', 'run_prefix',
                            'runid', 'sample_plate', 'sequencing_meth',
                            'linker', 'primer', 'primer_plate', 'well_id_384',
                            'plating', 'extractionkit_lot', 'extraction_robot',
                            'tm1000_8_tool', 'primer_date', 'mastermix_lot',
                            'water_lot', 'processing_robot', 'tm300_8_tool',
                            'tm50_8_tool', 'project_name', 'orig_name',
                            'well_description', 'pcr_primers', 'target_gene',
                            'tm10_8_tool', 'target_subfragment', 'well_id_96'}

    METAGENOMIC_PTYPE = 'Metagenomic'
    METATRANSCRIPTOMIC_PTYPE = 'Metatranscriptomic'
    AMPLICON_PTYPE = 'Amplicon'

    pipeline_types = {METAGENOMIC_PTYPE, AMPLICON_PTYPE,
                      METATRANSCRIPTOMIC_PTYPE}

    def __init__(self, configuration_file_path, run_id, sample_sheet_path,
                 mapping_file_path, output_path, qiita_job_id, pipeline_type):
        """
        Initialize Pipeline object w/configuration information.
        :param configuration_file_path: Path to configuration.json file.
        :param run_id: Used w/search_paths to locate input run_directory.
        :param sample_sheet_path: Path to sample-sheet.
        :param mapping_file_path: Path to mapping file.
        :param output_path: Path where all pipeline-generated files live.
        :param qiita_job_id: Qiita Job ID creating this Pipeline.
        :param pipeline_type: Pipeline type ('Amplicon', 'Metagenomic', etc.)
        """
        if sample_sheet_path is not None and mapping_file_path is not None:
            raise PipelineError("sample_sheet_path or mapping_file_path "
                                "must be defined, but not both.")

        if sample_sheet_path is None and mapping_file_path is None:
            raise PipelineError("sample_sheet_path or mapping_file_path "
                                "must be defined, but not both.")

        if pipeline_type not in Pipeline.pipeline_types:
            raise PipelineError(f"'{type}' is not a valid pipeline type.")

        self.pipeline_type = pipeline_type

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

        # If our extended validate() method discovers any warnings or
        # Errors, it will raise a PipelineError and return them w/in the
        # error message as a single string separated by '\n'.
        self.warnings = []
        self._directory_check(output_path, create=False)
        self.output_path = output_path
        self.run_id = run_id
        self.qiita_job_id = qiita_job_id
        self.pipeline = []

        if sample_sheet_path:
            self.search_paths = config['search_paths']
            self.sample_sheet = self._validate_sample_sheet(sample_sheet_path)
            self.mapping_file = None
        else:
            self.search_paths = config['amplicon_search_paths']
            self.mapping_file = self._validate_mapping_file(mapping_file_path)
            # unlike _validate_sample_sheet() which returns a SampleSheet
            # object that stores the path to the file it was created from,
            # _validate_mapping_file() just returns a DataFrame. Store the
            # path to the original mapping file itself as well.
            self.mapping_file_path = mapping_file_path
            self.sample_sheet = None

        self.run_dir = self._search_for_run_dir()

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

        if self.mapping_file is not None:
            # create dummy sample-sheet
            output_fp = join(output_path, 'dummy_sample_sheet.csv')
            self.generate_dummy_sample_sheet(self.run_dir, output_fp)
            self.sample_sheet = output_fp

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

    def run(self, callback=None):
        """
        Run all jobs added to Pipeline in the order they were added.
        :param callback: Optional function to call and upstate status with.
        :param callback(jid=): string identifying the current running process.
        :param callback(status=): a string message or description.
        :return:
        """
        for job in self.pipeline:
            job.run(callback=callback)

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
        :return: If successful, a valid sample-sheet. Raises descriptive
                 PipelineError() on all failures. Warning messages are
                 appended to self.warnings.
        """
        # validate the sample-sheet using metapool package.
        sheet = load_sample_sheet(sample_sheet_path)

        msgs = sheet.quiet_validate_and_scrub_sample_sheet()

        if any([isinstance(m, ErrorMessage) for m in msgs]):
            # msgs will contain both ErrorMessages and WarningMessages.
            # we want to identify if there are any messages and if so, create
            # a separate list for them. An Error should only be raised on
            # Error messages and in this case, all error messages should be
            # concatenated.
            errors = [x for x in msgs if isinstance(x, ErrorMessage)]

            if errors:
                msgs = [str(x).replace('ErrorMessage: ', '') for x in msgs]
                msgs = 'Sample-sheet contains errors:\n' + '\n'.join(msgs)
                raise PipelineError(msgs)
            else:
                raise PipelineError('Cannot parse sample-sheet.')
        else:
            # perform extended validation based on required fields for
            # seqpro, and other issues encountered.
            bioinformatics = sheet.Bioinformatics
            if 'library_construction_protocol' not in bioinformatics:
                msgs.append(ErrorMessage("column 'library_construction_protoco"
                                         "l' not found in Bioinformatics secti"
                                         "on"))
            if 'experiment_design_description' not in bioinformatics:
                msgs.append(ErrorMessage("column 'experiment_design_descriptio"
                                         "n' not found in Bioinformatics secti"
                                         "on"))

            if sheet.Header['Assay'] not in Pipeline.assay_types:
                msgs.append(ErrorMessage("Valid Assay values are "
                                         f"{Pipeline.assay_types}"))

            # look for duplicate samples. metapool will allow two rows w/the
            # same lane and sample_id if one or more other columns are
            # different. However seqpro expects the tuple (lane, sample_id) to
            # be unique for indexing.
            unique_indexes = []
            for item in sheet.samples:
                unique_index = f'{item.lane}_{item.sample_id}'
                if unique_index in unique_indexes:
                    msgs.append(ErrorMessage("A sample already exists with la"
                                             f"ne {item.lane} and sample-id "
                                             f"{item.sample_id}"))
                else:
                    unique_indexes.append(unique_index)

            errors = [x for x in msgs if isinstance(x, ErrorMessage)]

            if errors:
                msgs = [str(x).replace('ErrorMessage: ', '') for x in msgs]
                msgs = 'Sample-sheet contains errors:\n' + '\n'.join(msgs)
                raise PipelineError(msgs)

            # return a valid sample-sheet, and preserve any warning
            # messages
            self.warnings += [str(x) for x in msgs if
                              isinstance(x, WarningMessage)]
            return sheet

    def _validate_mapping_file(self, mapping_file_path):
        """
        Performs validation for mapping-files.
        :return: If successful, a valid mapping-file. Raises descriptive
                 PipelineError() on all failures. Warning messages are
                 appended to self.warnings.
        """
        try:
            df = pd.read_csv(mapping_file_path, delimiter='\t', dtype=str)
        except pd.errors.ParserError:
            raise PipelineError('Cannot parse mapping-file.')

        # first, detect any duplicate column names, regardless of any mixed-
        # capitalization, and notify the user.
        d = defaultdict(list)
        for column in df.columns:
            d[column.lower()].append(column)

        # generate a list of all unique column names that appear more than
        # once, regardless of capitalization. Then generate a list containing
        # lists of duplicate column names in their original case to report to
        # the user.
        dupes = [d[column] for column in
                 [col for col in d.keys() if len(d[col]) > 1]]

        if dupes:
            # column-names are case-insensitive, and must be unique.
            # return groups of duplicate column names (differentiated only by
            # a different mixed-case) to the user.
            raise PipelineError("Mapping-file contains duplicate columns: "
                                "%s" % ', '.join([str(tpl) for tpl in dupes]))

        # if columns are unique, determine if any columns are missing and/or
        # unexpected and notify the user.
        obs = set(df.columns.str.lower())

        # Note that Pipeline.mapping_file_columns is expected to be all lower-
        # case.

        # if an expected column is missing in observed, that is an error.
        # Note that since a mapping-file is just a DataFrame, there isn't a
        # distinction between a mapping-file that is missing n columns and has
        # n additional columns and a dataframe that is not a mapping-file at
        # all. This method assumes an external test has determined that the
        # file is a mapping-file already.
        missing_columns = Pipeline.mapping_file_columns - obs
        if missing_columns:
            raise PipelineError("Mapping-file is missing columns: "
                                "%s" % ', '.join(sorted(missing_columns)))

        # if an observed column is unexpected, that is a warning.
        unexpected_columns = obs - Pipeline.mapping_file_columns
        if unexpected_columns:
            self.warnings += [("Mapping-file contains additional columns: "
                               "%s" % ', '.join(unexpected_columns))]

        # rename all columns to their lower-case versions.
        # we will want to return this version to the user.
        df.columns = df.columns.str.lower()

        return df

    def generate_sample_info_files(self, addl_info=None):
        """
        Generate sample-information files in self.output_path.
        :param addl_info: A df of (sample-name, project-name) pairs.
        :return: A list of paths to sample-information-files.
        """
        if self.mapping_file is not None:
            # Generate a list of BLANKs for each project.
            df = self.mapping_file[['sample_name', 'project_name']]
        else:
            # Aggregate all data into a DataFrame
            data = [[x['Sample_ID'], x['Sample_Project']] for
                    x in self.sample_sheet.samples]
            df = pd.DataFrame(data, columns=['sample_name', 'project_name'])

        if addl_info is not None:
            df = pd.concat([df, addl_info],
                           ignore_index=True).drop_duplicates()

        df = df[df["sample_name"].str.startswith("BLANK") == True]  # noqa
        samples = list(df.to_records(index=False))
        projects = df.project_name.unique()

        paths = []
        for project in projects:
            samples_in_proj = [x for x, y in samples if y == project]
            some_path = join(self.output_path,
                             f'{self.run_id}_{project}_blanks.tsv')
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
                    row['title'] = sub(r'_\d+$', r'', project)

                    # generate values for the four columns that must be
                    # determined from sample-sheet information.

                    # convert 'BLANK14_10F' to 'BLANK14.10F', etc.
                    row['sample_name'] = sample.replace('_', '.')
                    row['host_subject_id'] = sample.replace('_', '.')
                    row['description'] = sample.replace('_', '.')
                    row['collection_timestamp'] = self.get_date_from_run_id()

                    row = [row[x] for x in Pipeline.sif_header]
                    f.write('\t'.join(row) + '\n')

        return paths

    def get_date_from_run_id(self):
        # assume all run_ids begin with coded datestamp:
        # 210518_...
        # allow exception if substrings cannot convert to int
        # or if array indexes are out of bounds.
        year = int(self.run_id[0:2]) + 2000
        month = int(self.run_id[2:4])
        day = int(self.run_id[4:6])
        return f'{year}-{month}-{day}'

    def get_sample_ids(self):
        '''
        Returns list of sample-ids sourced from sample-sheet or pre-prep file
        :return: list of sample-ids
        '''

        # test for self.mapping_file, since self.sample_sheet will be
        # defined in both cases.
        if self.mapping_file is not None:
            results = list(self.mapping_file.sample_name)
        else:
            results = [x.Sample_ID for x in self.sample_sheet.samples]

        return results

    def get_sample_names(self, project_name=None):
        '''
        Returns list of sample-names sourced from sample-sheet or pre-prep file
        :param project_name: If None, return all sample-names.
        :return: list of sample-names
        '''
        # test for self.mapping_file, since self.sample_sheet will be
        # defined in both cases.
        if self.mapping_file is not None:
            return self._get_sample_names_from_mapping_file(project_name)
        else:
            return self._get_sample_names_from_sample_sheet(project_name)

    def _get_sample_names_from_sample_sheet(self, project_name):
        if project_name is None:
            return [x.Sample_Name for x in self.sample_sheet.samples]
        else:
            # Since the project-name is stored in an internal variable
            # in a third-party library, convert the data structure to
            # JSON using the exposed method and obtain from the result.
            jsn = json_loads(self.sample_sheet.to_json())
            return [x['Sample_Name'] for x in jsn['Data'] if
                    f'{project_name}_' in x['Sample_Project']]

    def _get_sample_names_from_mapping_file(self, project_name):
        if project_name is None:
            return list(self.mapping_file.sample_name)
        else:
            df = self.mapping_file[self.mapping_file['project_name'] ==
                                   project_name]
            return list(df['sample_name'])

    def _parse_project_name(self, project_name, short_names):
        '''
        Split fully-qualified project_name into a project_name and a qiita-id
        if possible. Else return project_name and None.
        :param project_name: A fully-qualified project name e.g: Feist_1161.
        :param short_names: True returns orig. value. False returns name only.
        :return: Tuple (project-name, qiita-id)
        '''
        # This functionality could be folded into metapool package's
        # remove_qiita_id() in the future.
        if project_name is None or project_name == '':
            raise ValueError("project_name cannot be None or empty string")

        matches = search(r'^(.+)_(\d+)$', str(project_name))

        if matches is None:
            raise ValueError(f"'{project_name}' does not contain a Qiita-ID.")

        if short_names is False:
            # return the fully-qualified project name w/Qiita ID.
            return project_name, matches[2]
        else:
            # return the project's name and qiita_id
            return matches[1], matches[2]

    def get_project_info(self, short_names=False):
        # test for self.mapping_file, since self.sample_sheet will be
        # defined in both cases.
        results = []

        if self.mapping_file is not None:
            sample_project_map = {pn: _df.sample_name.values for pn, _df in
                                  self.mapping_file.groupby('project_name')}

            for project in sample_project_map:
                p_name, q_id = self._parse_project_name(project, short_names)
                results.append(
                    {'project_name': p_name, 'qiita_id': q_id})
        else:
            bioinformatics = self.sample_sheet.Bioinformatics
            for res in bioinformatics.to_dict('records'):
                p_name, q_id = self._parse_project_name(res['Sample_Project'],
                                                        short_names)
                results.append({'project_name': p_name,
                                'qiita_id': q_id})

        return results

    @staticmethod
    def is_mapping_file(mapping_file_path):
        '''
        Returns True if file follows basic mapping-file format.
        '''
        try:
            df = pd.read_csv(mapping_file_path, delimiter='\t', dtype=str)
        except pd.errors.ParserError:
            return False

        # if the expected subset of columns required for a mapping-file
        # are present, then consider this a mapping file, even if it's
        # an invalid one.
        exp_columns = frozenset({'barcode', 'tm1000_8_tool',
                                 'extraction_robot', 'pcr_primers'})

        return set(df.columns.str.lower()).issuperset(exp_columns)

    @staticmethod
    def is_sample_sheet(sample_sheet_path):
        '''
        Returns True if file follows basic sample-sheet format.
        '''

        # Check to see if the file begins w/[Header].
        # Ignoring any legacy comments.
        with open(sample_sheet_path, 'r') as f:
            line = f.readline()
            while line:
                if line.startswith('#'):
                    line = f.readline()
                else:
                    break

            if line.startswith('[Header]'):
                return True

        return False

    def _generate_dummy_sample_sheet(self, first_read, last_read,
                                     indexed_reads, dummy_sample_id):
        # create object and initialize header
        sheet = AmpliconSampleSheet()
        sheet.Header['IEMFileVersion'] = '4'
        sheet.Header['Date'] = '10/27/22'
        sheet.Header['Workflow'] = 'GenerateFASTQ'
        sheet.Header['Application'] = 'FASTQ Only'
        sheet.Header['Assay'] = 'TruSeq HT'
        sheet.Header['Description'] = 'test_run'
        sheet.Header['Chemistry'] = 'Amplicon'

        # generate override_cycles string
        tmp = [f"N{x['NumCycles']}" for x in indexed_reads]
        tmp = ';'.join(tmp)
        override_cycles = f"Y{first_read};{tmp};Y{last_read}"

        # set Reads and Settings according to input values
        # we'll get this from the code on the server
        sheet.Reads = [first_read, last_read]
        sheet.Settings['OverrideCycles'] = override_cycles
        sheet.Settings['MaskShortReads'] = '1'
        sheet.Settings['CreateFastqForIndexReads'] = '1'

        dummy_samples = {'Sample_ID': dummy_sample_id,
                         'Sample_Plate': '',
                         'Sample_Well': '',
                         'I7_Index_ID': '',
                         'index': '',
                         'I5_Index_ID': '',
                         'index2': ''
                         }
        sheet.add_sample(sample_sheet.Sample(dummy_samples))

        # contacts won't matter for the dummy sample-sheet.
        contacts = [['c2cowart@ucsd.edu', 'SomeProject'],
                    ['antgonza@gmail.com', 'AnotherProject']]

        # we'll get these from input parameters as well.
        contacts = pd.DataFrame(columns=['Email', 'Sample_Project'],
                                data=contacts)
        sheet.Contact = contacts

        # add a dummy sample.
        samples = [[dummy_sample_id, 'NA', 'NA',
                    'FALSE', 'FALSE', '14782']]

        samples = pd.DataFrame(columns=['Project', 'ForwardAdapter',
                                        'ReverseAdapter', 'PolyGTrimming',
                                        'HumanFiltering', 'QiitaID'],
                               data=samples)

        sheet.Bioinformatics = samples

        return sheet

    def generate_dummy_sample_sheet(self, run_dir, output_fp):
        if exists(run_dir):
            reads = self.process_run_info_file(join(run_dir, 'RunInfo.xml'))
        else:
            raise ValueError("run_dir %s not found." % run_dir)

        # assumptions are first and last reads are non-indexed and there
        # are always two. Between them there is either 1 or 2 indexed
        # reads. If this is not true, raise an Error.

        if len(reads) < 3 or len(reads) > 4:
            # there must be a first and last read w/a minimum of one read
            # in the middle and maximum two in the middle.
            raise ValueError("RunInfo.xml contains abnormal reads.")

        first_read = reads.pop(0)
        last_read = reads.pop()

        if (first_read['IsIndexedRead'] is True or
                last_read['IsIndexedRead'] is True):
            raise ValueError("RunInfo.xml contains abnormal reads.")

        # confirm the interior read(s) are indexed ones.
        for read in reads:
            if read['IsIndexedRead'] is False:
                raise ValueError("RunInfo.xml contains abnormal reads.")

        dummy_sample_id = basename(run_dir) + '_SMPL1'

        sheet = self._generate_dummy_sample_sheet(first_read['NumCycles'],
                                                  last_read['NumCycles'],
                                                  reads, dummy_sample_id)

        with open(output_fp, 'w') as f:
            sheet.write(f, 1)

    def process_run_info_file(self, run_info_fp):
        def process_reads(reads):
            # extract all read elements as a list.
            # the contents of each Read element are highly regular.
            # for now, process w/out installing xml2dict or other
            # library into Qiita env.
            found = findall('<Read (.+?) />', reads)

            results = []
            for item in found:
                attributes = item.split(' ')
                d = {}
                for attribute in attributes:
                    k, v = attribute.split('=')
                    if k in ['NumCycles', 'Number']:
                        v = int(v.strip('"'))
                    elif k in ['IsIndexedRead']:
                        v = v.strip('"')
                        v = False if v == 'N' else True
                    else:
                        raise ValueError("Unknown key: %s" % k)
                    d[k] = v
                results.append(d)

            return results

        with open(run_info_fp, 'r') as f:
            s = f.read()
            reads = search('<Reads>(.+?)</Reads>', s.replace('\n', ''))
            if reads:
                result = reads.group(1)
            else:
                raise ValueError("Cannot extract read information")
            return process_reads(result)
