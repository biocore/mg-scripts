from json import load as json_load
from json.decoder import JSONDecodeError
from os import makedirs, listdir
from os.path import join, exists, isdir, basename
from metapool import KLSampleSheet, quiet_validate_and_scrub_sample_sheet
from metapool.plate import ErrorMessage, WarningMessage
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
import logging
from re import sub, findall, search
import sample_sheet
import pandas as pd


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
                 mapping_file_path, output_path, qiita_job_id,
                 config_dict=None):
        """
        Initialize Pipeline object w/configuration information.
        :param configuration_file_path: Path to configuration.json file.
        :param run_id: Used w/search_paths to locate input run_directory.
        :param sample_sheet_path: Path to sample-sheet.
        :param mapping_file_path: Path to mapping file.
        :param output_path: Path where all pipeline-generated files live.
        :param qiita_job_id: Qiita Job ID creating this Pipeline.
        :param config_dict: (Optional) Dict used instead of config file.
        """
        if sample_sheet_path is not None and mapping_file_path is not None:
            raise PipelineError("sample_sheet_path or mapping_file_path "
                                "must be defined, but not both.")

        if sample_sheet_path is None and mapping_file_path is None:
            raise PipelineError("sample_sheet_path or mapping_file_path "
                                "must be defined, but not both.")

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
        :param callback(id=): a string identifying the current running process.
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
        if self.mapping_file is not None:
            df = self.mapping_file[['sample_name', 'project_name']]
            samples = list(df.to_records(index=False))
            samples = [x for x in samples if x[0].startswith('BLANK')]
            projects = list(set([y for x, y in samples]))
        else:
            samples = []
            for sample in self.sample_sheet.samples:
                if sample['Sample_ID'].startswith('BLANK'):
                    samples.append((sample['Sample_ID'],
                                    sample['Sample_Project']))
            projects = list(set([y for x, y in samples]))

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
        # test for self.mapping_file, since self.sample_sheet will be
        # defined in both cases.
        if self.mapping_file is not None:
            results = []
            sample_project_map = {pn: _df.sample_name.values for pn, _df in
                                  self.mapping_file.groupby('project_name')}

            for project in sample_project_map:
                results += sample_project_map[project]
        else:
            results = [x.Sample_ID for x in self.sample_sheet.samples]

        return results

    def get_project_info(self):
        # test for self.mapping_file, since self.sample_sheet will be
        # defined in both cases.
        results = []

        if self.mapping_file is not None:
            sample_project_map = {pn: _df.sample_name.values for pn, _df in
                                  self.mapping_file.groupby('project_name')}

            for project in sample_project_map:
                qiita_id = project.split('_')[-1]
                results.append(
                    # assume project names end in qiita ids
                    {'project_name': project, 'qiita_id': qiita_id})
        else:
            bioinformatics = self.sample_sheet.Bioinformatics
            for result in bioinformatics.to_dict('records'):
                results.append({'project_name': result['Sample_Project'],
                                'qiita_id': result['QiitaID']})

        return results

    @staticmethod
    def is_mapping_file(mapping_file_path):
        try:
            Pipeline._validate_mapping_file(mapping_file_path)
            return True
        except PipelineError:
            return False

    @staticmethod
    def _validate_mapping_file(mapping_file_path):
        exp = {'barcode', 'library_construction_protocol', 'mastermix_lot',
               'sample_plate', 'center_project_name', 'instrument_model',
               'tm1000_8_tool', 'well_id', 'tm50_8_tool', 'well_description',
               'run_prefix', 'run_date', 'center_name', 'tm300_8_tool',
               'extraction_robot', 'experiment_design_description',
               'platform', 'water_lot', 'project_name', 'pcr_primers',
               'sequencing_meth', 'plating', 'orig_name', 'linker', 'runid',
               'target_subfragment', 'primer', 'primer_plate', 'sample_name',
               'run_center', 'primer_date', 'target_gene', 'processing_robot',
               'extractionkit_lot'}
        try:
            df = pd.read_csv(mapping_file_path, delimiter='\t')

            # if the two sets of headers are equal, then return the dataframe
            # and no error message.
            if len(exp - set(df.columns)) == 0:
                return df

            msg = 'missing columns: %s' % ', '.join(exp - set(df.columns))
        except pd.errors.ParserError:
            # ignore parser errors as they obviously prove this is not a
            # valid mapping file.
            msg = 'could not parse file "%s"' % mapping_file_path

        raise PipelineError(msg)

    def _generate_dummy_sample_sheet(self, first_read, last_read,
                                     indexed_reads, dummy_sample_id):
        # create object and initialize header
        sheet = KLSampleSheet()
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

    def get_sample_project_map(self, mapping_file_df):
        sample_project_map = {pn: _df.sample_name.values for pn, _df in
                              self.mapping_file.groupby('project_name')}

        for sample_name, project_name in zip(mapping_file_df.sample_name,
                                             mapping_file_df.project_name):
            if project_name not in sample_project_map:
                sample_project_map[project_name] = []
            sample_project_map[project_name].append(sample_name)

        return sample_project_map
