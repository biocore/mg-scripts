from json import load as json_load
from json.decoder import JSONDecodeError
from os import makedirs, listdir
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


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)


class Pipeline:
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
        self.sentinel_file = "RTAComplete.txt"
        self.qiita_job_id = qiita_job_id
        self.pipeline = []

        if self.older_than < 0 or self.younger_than < 0:
            raise PipelineError('older_than and younger_than cannot be '
                                'less than zero.')

        if self.older_than >= self.younger_than:
            raise PipelineError(
                'older_than cannot be equal to or less than younger_than.')

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


if __name__ == '__main__':
    logging.debug("Starting Pipeline")

    sample_sheet_path = '/path/to/sample-sheet'
    run_id = '210518_A00953_0305_test11'
    config_file_path = '/path/to/configuration.json'
    qiita_job_id = 'NOT_A_QIITA_JOB_ID'
    output_path = '/path/to/output_path'
    config_dict = None

    pipeline = Pipeline(config_file_path, run_id, output_path, config_dict,
                        qiita_job_id)

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

    pipeline.add(QCJob(pipeline.run_dir,
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
                   config['job_pool_size']))

    config = pipeline.configuration['fastqc']

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
                           config['multiqc_config_file_path']))

    config = pipeline.configuration['seqpro']

    pipeline.add(GenPrepFileJob(pipeline.run_dir,
                                pipeline.output_path,
                                sdo.sample_sheet_path,
                                config['seqpro_path'],
                                config['modules_to_load'],
                                qiita_job_id))

    pipeline.run()

    # Currently this is performed outside a Job as it's not the
    # responsibility of ConvertJob() to see that this file is copied for
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

    # TODO: copy the results
