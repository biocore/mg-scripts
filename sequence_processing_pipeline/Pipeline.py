from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.QCJob import QCJob
# from sequence_processing_pipeline.FastQC import FastQCJOb
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.PipelineError import PipelineError
from time import time as epoch_time
import logging
import os
from os.path import join, exists
import json
from json.decoder import JSONDecodeError


class Pipeline:
    def __init__(self, configuration_file_path, input_directory,
                 output_directory):

        try:
            f = open(configuration_file_path)
            self.configuration = json.load(f)
        except FileNotFoundError:
            raise PipelineError(f'{configuration_file_path} does not exist.')
        except JSONDecodeError:
            raise PipelineError(f'{configuration_file_path} is not a valid ',
                                'json file')

        config = self.configuration['pipeline']
        younger_than = config['younger_than']
        older_than = config['older_than']
        archive_path = config['archive_path']

        if output_directory == archive_path:
            raise PipelineError(
                "output_directory '%s' is the same as "
                "archive_path '%s'." % (
                    output_directory, archive_path))

        self._directory_check(input_directory, create=False)
        self._directory_check(output_directory, create=True)
        self._directory_check(archive_path, create=True)

        '''
        if nprocs > 16:
            raise PipelineError('nprocs cannot exceed 16.')
        elif nprocs < 1:
            raise PipelineError('nprocs cannot be less than 1.')
        '''

        if older_than >= younger_than:
            raise PipelineError(
                'older_than cannot be equal to or less than younger_than.')

        self.run_dir = input_directory

        logging.debug("Root directory: %s" % self.run_dir)

        self.sentinel_file = "RTAComplete.txt"
        logging.debug("Sentinel File Name: %s" % self.sentinel_file)
        # internally, threshold will be represented in seconds
        self.younger_than = younger_than * 60 * 60
        self.older_than = older_than * 60 * 60
        s = "Filter directories younger than %d hours old" % younger_than
        logging.debug(s)
        s = "Filter directories older than %d hours old" % older_than
        logging.debug(s)
        # self.nprocs = nprocs
        # logging.debug("nprocs: %d" % nprocs)
        self.output_dir = output_directory
        logging.debug("Output Directory: %s" % self.output_dir)
        # location in /sequencing directory or other directory
        # must end in name of the form 200nnn_xnnnnn_nnnn_xxxxxxxxxx or
        # similar.
        self.final_output_dir = archive_path
        logging.debug("Final Output Directory: %s" % self.final_output_dir)

    def _directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    os.makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PipelineError, so it gets handled in the same location
                    # as the others.
                    raise PipelineError(str(e))
            else:
                raise PipelineError("directory_path '%s' does not exist." %
                                    directory_path)

    def process(self, sample_sheet_path):
        '''
        Process a single BCL directory.
        Assume sample sheet is supplied externally.
        Assume directory doesn't require time period checking.
        :param sample_sheet:
        :return:
        '''
        try:
            sdo = SequenceDirectory(self.run_dir,
                                    external_sample_sheet=sample_sheet_path)
            sample_sheet_params = sdo.process()

            ss_path = sample_sheet_params['sample_sheet_path']
            fastq_output_directory = join(self.run_dir, 'Data', 'Fastq')

            # TODO: Add checks to ensure fields are defined.
            config = self.configuration['bcl-convert']
            convert_job = ConvertJob(self.run_dir,
                                     ss_path,
                                     fastq_output_directory,
                                     True,
                                     config['queue'],
                                     config['nodes'],
                                     config['nprocs'],
                                     config['wallclock_time_in_hours'] * 3600,
                                     config['per_process_memory_limit'])

            convert_job.run()

            config = self.configuration['qc']
            qc_job = QCJob(self.run_dir,
                           sample_sheet_params['sample_sheet_path'],
                           config['mmi_db'],
                           self.final_output_dir,
                           config['queue'],
                           config['nodes'],
                           config['nprocs'],
                           config['wallclock_time_in_hours'] * 3600,
                           config['per_process_memory_limit'],
                           config['fastp_executable_path'],
                           config['minimap2_executable_path'],
                           config['samtools_executable_path'])

            qc_job.run()

        except PipelineError as e:
            logging.error(e)

    def _time_is_right(self, timestamp):
        # calculate how old the timestamp is
        delta_t = epoch_time() - timestamp

        # if the timestamp falls w/in the range defined by the object, then
        # it is a legitimate target.
        if delta_t < self.younger_than and delta_t > self.older_than:
            return True

        return False

    def _find_bcl_directories(self):
        """
        Walk root directory and locate all scan directory root folders
        :return: list of BCL directories within root folder.
        """
        new_dirs = []

        for root, dirs, files in os.walk(self.run_dir):
            for some_directory in dirs:
                some_path = join(root, some_directory)
                if '/Data/Intensities/BaseCalls/' in some_path:
                    # the root directory of every scan directory will have
                    # this substring in the path of one or more of their
                    # subdirectories. By collecting these subdirectories
                    # and extracting the root directory of each one, we can
                    # build a set of unique root directories.
                    s = some_path.split('/Data/Intensities/BaseCalls')[0]
                    some_path = s
                    new_dirs.append(some_path)

        # remove duplicates
        new_dirs = list(set(new_dirs))

        for new_dir in new_dirs:
            logging.info("%s found." % new_dir)

        return new_dirs

    def _filter_directories_for_time(self, new_dirs):
        """
        Filter directories for those that match allowed timespan.
        :return: list of BCL directories within timespan.
        """
        filtered_dirs = []
        for some_path in new_dirs:
            # whether a single BCL directory is passed, or a nested tree
            # of directories is passed, assume that a potential BCL
            # directory must contain a file named self.sentinel_file.
            if os.path.exists(join(some_path, self.sentinel_file)):
                some_timestamp = os.path.getmtime(some_path)
                # if the last modified timestamp of the directory is
                # within the threshold of what is considered 'new and
                # completed data',then this directory is a legitimate
                # target.
                if self._time_is_right(some_timestamp):
                    # save the path as well as the original epoch timestamp
                    # as a tuple.
                    filtered_dirs.append((some_path, some_timestamp))
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
