from json import load as json_load
from json.decoder import JSONDecodeError
from os import makedirs, listdir, walk
from os.path import join, exists, isdir, getmtime
from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from time import time as epoch_time
import logging
import shutil

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


class Pipeline:
    def __init__(self, configuration_file_path, run_id, config_dict=None):
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
                                    "{self.configuration_file_path}")

        self.search_paths = config['search_paths']
        self.run_id = run_id
        self.run_dir = self._search_for_run_dir()
        self.products_dir = join(self.run_dir, run_id)
        makedirs(self.products_dir, exist_ok=True)
        self.fastq_output_dir = join(self.run_dir, 'Data', 'Fastq')
        self.final_output_dir = config['archive_path']
        self.younger_than = config['younger_than']
        self.older_than = config['older_than']
        if self.older_than < 0 or self.younger_than < 0:
            raise PipelineError('older_than and younger_than cannot be '
                                'less than zero.')
        if self.older_than >= self.younger_than:
            raise PipelineError(
                'older_than cannot be equal to or less than younger_than.')
        # self.final_output_dir needs to already exist.
        # join(self.final_output_dir, self.run_id) should not already exist.
        self.directory_check(self.final_output_dir, create=False)
        tmp_path = join(self.final_output_dir, self.run_id)
        if exists(tmp_path):
            raise PipelineError(f"{tmp_path} already exists")
        self.run_id = run_id
        self.sentinel_file = "RTAComplete.txt"

    def _search_for_run_dir(self):
        for search_path in self.search_paths:
            logging.debug(f'Searching {search_path} for {self.run_id}')
            for entry in listdir(search_path):
                some_path = join(search_path, entry)
                # ensure some_path never ends in '/'
                some_path = some_path.rstrip('/')
                if isdir(some_path) and some_path.endswith(self.run_id):
                    logging.debug(f'Found {some_path}')
                    return some_path
        raise PipelineError(f"A run-dir for '{self.run_id}' could not be "
                            "found")

    def copy_results_to_archive(self):
        # hack an empty Job() object to use job._system_call().
        job = Job('', '', [])
        job._system_call('rsync -avp %s %s' % (self.products_dir,
                                               self.final_output_dir))

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

    def find_bcl_directories(self, root_path):
        """
        Walk root directory and locate all scan directory root folders
        :return: list of BCL directories within root folder.
        """
        new_dirs = []

        for root, dirs, files in walk(root_path):
            for some_directory in dirs:
                some_path = join(root, some_directory)
                if '/Data/Intensities/BaseCalls/' in some_path:
                    # the root directory of every scan directory will have
                    # this substring in the path of one or more of their
                    # subdirectories. By collecting these subdirectories
                    # and extracting the root directory of each one, we can
                    # build a set of unique root directories.
                    new_dirs.append(
                        some_path.split('/Data/Intensities/BaseCalls')[0])

        # remove duplicates
        new_dirs = list(set(new_dirs))

        return new_dirs

    def filter_directories_for_time(self, new_dirs):
        """
        Filter directories for those that match allowed timespan.
        :return: list of BCL directories within timespan.
        """
        filtered_dirs = []
        for some_path in new_dirs:
            # whether a single BCL directory is passed, or a nested tree
            # of directories is passed, assume that a potential BCL
            # directory must contain a file named self.sentinel_file.
            if exists(join(some_path, self.sentinel_file)):
                some_timestamp = getmtime(some_path)
                # if the last modified timestamp of the directory is
                # within the threshold of what is considered 'new and
                # completed data',then this directory is a legitimate
                # target.
                delta_t = epoch_time() - some_timestamp
                if delta_t < (self.younger_than * 3600) and delta_t > (
                        self.older_than * 3600):
                    # save the path as well as the original epoch timestamp
                    # as a tuple.
                    filtered_dirs.append(some_path)
                else:
                    s = "The timestamp for {} is not within bounds."
                    logging.debug(s.format(some_path))
            else:
                # This is a warning, rather than an error because a
                # directory of BCL directories would be a valid parameter,
                # even though it doesn't contain BCL data itself.
                s = "{} does not contain a file named '{}'."
                logging.warning(s.format(some_path, self.sentinel_file))

        return filtered_dirs


if __name__ == '__main__':
    logging.debug("Starting Pipeline")

    sample_sheet_path = '/path/to/sample-sheet'
    run_id = '210518_A00953_0305_test11'
    config_file_path = '/path/to/configuration.json'
    qiita_job_id = 'NOT_A_QIITA_JOB_ID'

    pipeline = Pipeline(config_file_path, run_id)
    sdo = SequenceDirectory(pipeline.run_dir,
                            sample_sheet_path=sample_sheet_path)

    config = pipeline.configuration['bcl-convert']
    convert_job = ConvertJob(pipeline.run_dir,
                             sdo.sample_sheet_path,
                             pipeline.fastq_output_dir,
                             config['queue'],
                             config['nodes'],
                             config['nprocs'],
                             config['wallclock_time_in_hours'],
                             config['per_process_memory_limit'],
                             config['executable_path'],
                             config['modules_to_load'],
                             qiita_job_id)

    logging.debug("Starting ConvertJob")
    convert_job.run()
    logging.debug("ConvertJob Finished")

    # Currently this is performed outside a Job as it's not the
    # resposibility of ConvertJob() to see that this file is copied for
    # GenPrepFileJob(), and GenPrepFileJob() shouldn't have to know where
    # this information is located.
    src1 = join(pipeline.run_dir, 'Data/Fastq/Reports')
    src2 = join(pipeline.run_dir, 'Data/Fastq/Stats')
    if exists(src1):
        shutil.copytree(src1, join(pipeline.products_dir, 'Reports'))
    elif exists(src2):
        shutil.copytree(src2, join(pipeline.products_dir, 'Stats'))
    else:
        raise PipelineError("Cannot locate Fastq metadata directory.")

    config = pipeline.configuration['qc']
    qc_job = QCJob(pipeline.run_dir,
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
                   pipeline.products_dir, )

    logging.debug("Starting QCJob")
    qc_job.run()
    logging.debug("QCJob Finished")

    config = pipeline.configuration['fastqc']
    output_directory = join(pipeline.products_dir, 'FastQC')
    fastqc_job = FastQCJob(pipeline.run_dir,
                           output_directory,
                           config['nprocs'],
                           config['nthreads'],
                           config['fastqc_executable_path'],
                           config['modules_to_load'],
                           qiita_job_id,
                           run_id,
                           config['queue'],
                           config['nodes'],
                           config['wallclock_time_in_hours'],
                           config['job_total_memory_limit'],
                           config['job_pool_size'])

    logging.debug("Starting FastQCJob")
    fastqc_job.run()
    logging.debug("FastQCJob Finished")

    gpf_output_path = join(pipeline.products_dir, 'prep-files')
    makedirs(gpf_output_path, exist_ok=True)

    config = pipeline.configuration['seqpro']
    gpf_job = GenPrepFileJob(pipeline.products_dir,
                             sdo.sample_sheet_path,
                             gpf_output_path,
                             config['seqpro_path'],
                             config['modules_to_load'],
                             qiita_job_id)

    try:
        logging.debug("Starting GenPrepFileJob")
        gpf_job.run()
        logging.debug("GenPrepFileJob Finished")
    except PipelineError as e:
        print(f"Caught known seqpro error: {str(e)}")

    logging.debug("Starting rsync")
    pipeline.copy_results_to_archive()
    logging.debug("Rsync Finished")
