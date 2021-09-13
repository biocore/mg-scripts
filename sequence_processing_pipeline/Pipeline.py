from sequence_processing_pipeline.BCLConvertJob import BCLConvertJob
from sequence_processing_pipeline.HumanFilterJob import HumanFilterJob
from sequence_processing_pipeline.FastQC import FastQCJOb
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.PipelineError import PipelineError
from time import time as epoch_time
import logging
import os


class Pipeline:
    # TODO: We need to change final_output_directory so that it if there are
    #  multiple output directories they can all be under output_directory.
    #  Really these are the outputs of stages.
    def __init__(self, input_directory, output_directory,
                 final_output_directory, younger_than=48,
                 older_than=24, nprocs=16, should_filter=False):
        self.root_dir = input_directory
        logging.debug("Root directory: %s" % self.root_dir)
        self.sentinel_file = "RTAComplete.txt"
        logging.debug("Sentinel File Name: %s" % self.sentinel_file)
        # internally, threshold will be represented in seconds
        self.younger_than = younger_than * 60 * 60
        self.older_than = older_than * 60 * 60
        s = "Filter directories younger than %d hours old" % younger_than
        logging.debug(s)
        s = "Filter directories older than %d hours old" % older_than
        logging.debug(s)
        self.should_filter = should_filter
        self.nprocs = nprocs
        self.output_dir = output_directory
        self.final_output_dir = final_output_directory

        # hard-code for now
        self.job_owner_home = '/home/jdereus'
        self.email_list = ['jdereus@health.ucsd.edu', 'ccowart@ucsd.edu']
        self.bclconvert_template = "{}/seq_proc_dev/bclconvert_human_slurm.sh".format(self.job_owner_home)

    def process(self):
        '''
        Process a path containing multiple BCL directories.
        Assume that sample sheets are stored w/in each directory.
        Assume that filtering directories using timestamps will detect
         the new directories.
        :return:
        '''
        directories_found = self._find_bcl_directories()
        results = self._filter_directories_for_time(directories_found)
        filtered_dirs_w_timestamps = results

        for sequence_directory, timestamp in filtered_dirs_w_timestamps:
            try:
                sdo = SequenceDirectory(sequence_directory)
                sample_sheet_params = sdo.process()
                s = sample_sheet_params['sample_sheet_path']

                self._generate_run_config_file(s)

                bcl_convert_job = BCLConvertJob(self.root_dir,
                                                self.should_filter)

                bcl_convert_job.run(sample_sheet_params['sequence_directory'],
                                    sample_sheet_params['sample_sheet_path'],
                                    sample_sheet_params['base_mask'],
                                    sample_sheet_params['experiment_name'],
                                    self.bclconvert_template)

                chemistry = sample_sheet_params['chemistry']

                human_filter_job = HumanFilterJob(sdo,
                                                  self.nprocs,
                                                  self.job_owner_home,
                                                  self.email_list,
                                                  chemistry,
                                                  self.output_dir,
                                                  self.final_output_dir)

                p1 = sample_sheet_params['sample_sheet_path']
                p2 = sample_sheet_params['sequence_directory']
                human_filter_job.run(p1, p2)

                output_dir = 'foo'
                project = 'foo'

                fast_qc_job = FastQCJOb(self.root_dir, output_dir, self.nprocs, project)

                fast_qc_job.run()

            except PipelineError as e:
                logging.error(e)

    def process_one(self, sample_sheet_path):
        '''
        Process a single BCL directory.
        Assume sample sheet is supplied externally.
        Assume directory doesn't require time period checking.
        :param sample_sheet:
        :return:
        '''
        try:
            sdo = SequenceDirectory(self.root_dir,
                                    external_sample_sheet=sample_sheet_path)
            sample_sheet_params = sdo.process()
            s = sample_sheet_params['sample_sheet_path']

            self._generate_run_config_file(s)

            bcl_convert_job = BCLConvertJob(self.root_dir,
                                            self.should_filter)

            bcl_convert_job.run(sample_sheet_params['sequence_directory'],
                                sample_sheet_params['sample_sheet_path'],
                                sample_sheet_params['base_mask'],
                                sample_sheet_params['experiment_name'],
                                self.bclconvert_template)

            chemistry = sample_sheet_params['chemistry']

            human_filter_job = HumanFilterJob(sdo,
                                              self.nprocs,
                                              self.job_owner_home,
                                              self.email_list,
                                              chemistry,
                                              self.output_dir,
                                              self.final_output_dir)

            p1 = sample_sheet_params['sample_sheet_path']
            p2 = sample_sheet_params['sequence_directory']
            human_filter_job.run(p1, p2)

            output_dir = 'foo'
            project = 'foo'

            fast_qc_job = FastQCJOb(self.root_dir, output_dir, self.nprocs, project)

            fast_qc_job.run()

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

        for root, dirs, files in os.walk(self.root_dir):
            for some_directory in dirs:
                some_path = join_path(root, some_directory)
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
            if os.path.exists(join_path(some_path, self.sentinel_file)):
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

    def _generate_run_config_file(self, sample_sheet_path):
        '''
        Generates a run_config.txt file in self.root based on data found in
        sample_sheet_path.
        :param sample_sheet_path: Path to a sample sheet CSV file.
        :return: None
        '''
        with open(sample_sheet_path, 'r') as f:
            # we can (and should) open up this file as a proper INI
            # file. However it may have trailing ',' characters,
            # which may impede this. Hence, open it as a regular
            # file and clean it up before processing.
            lines = f.readlines()
            # first, let's strip out \n and \r\n and any leading or
            # trailing whitespaces
            lines = [x.strip() for x in lines]
            # second, let's strip out trailing ',' characters, if
            # they're present. We'll keep any leading ',' characters
            # or in-line ',' characters.
            lines = [x.rstrip(',') for x in lines]
            # lastly, there have been cases where lines contain only
            # ',' characters. Hence, look for and remove these now
            # empty lines before continuing.
            lines = [x for x in lines if x]
            # since the file is already in memory, we won't use INI
            # library to parse it, (for now). Although it would be
            # cleaner.
            metadata = []
            sentinel = False
            for i in range(0, len(lines)):
                if lines[i] == '[Bioinformatics]':
                    # when Bioinformatics section is found, start
                    # copying lines to the list buffer, but don't
                    # copy the header itself.
                    sentinel = True
                elif lines[i].startswith('['):
                    # when the next header is found, stop copying
                    # lines to the list buffer. Don't include this
                    # header line, either.
                    sentinel = False
                elif sentinel is True:
                    # this should be a line in between
                    # [Bioinformatics] and the next section. Copy it
                    # to the buffer.
                    metadata.append(lines[i])
                # if the end of the file is reached before the next
                # header is found, that means there were no more
                # headers and that's okay too.

            # remove duplicate lines (this appears to be an issue in
            # the original bash scripts.)
            metadata = list(set(metadata))
            metadata.sort()

            # write the sanitized data out to legacy file.
            s = join_path(self.root_dir, 'run_config.txt')
            run_config_file_path = s
            with open(run_config_file_path, 'w') as f2:
                for line in metadata:
                    # end lines w/proper UNIX-style newline, unless
                    # Windows '\r\n' is preferred.
                    f2.write("%s\n" % line)

    def _get_contact_information(self, sample_sheet_path):
        '''
        Extracts contact information from sample sheet.
        :param sample_sheet_path:
        :return:
        '''
        with open(sample_sheet_path, 'r') as f:
            lines = f.readlines()
            lines = [x.strip() for x in lines]
            lines = [x.rstrip(',') for x in lines]
            lines = [x for x in lines if x]
            metadata = []
            sentinel = False
            for i in range(0, len(lines)):
                if lines[i] == '[Contact]':
                    sentinel = True
                elif lines[i].startswith('['):
                    sentinel = False
                elif sentinel is True:
                    metadata.append(lines[i])

            metadata = list(set(metadata))
            metadata.sort()

            return metadata

    def _get_column_count(self, sample_sheet_path):
        '''
        Extracts number of columns from sample sheet.
        This version returns the number of values, rather than the number of
        delimiters.
        :param sample_sheet_path:
        :return:
        '''
        with open(sample_sheet_path, 'r') as f:
            lines = f.readlines()
            lines = [x.strip() for x in lines]
            lines = [x.rstrip(',') for x in lines]
            lines = [x for x in lines if x]
            sentinel = False
            for i in range(0, len(lines)):
                if lines[i] == '[Data]':
                    sentinel = True
                elif lines[i].startswith('['):
                    sentinel = False
                elif sentinel is True:
                    count = len(lines[i].split(','))
                    logging.debug("Column Count (+1 over legacy): %s" % count)
                    return count

    def _get_reads(self, sample_sheet_path):
        '''
        Extracts contact information from sample sheet.
        :param sample_sheet_path:
        :return:
        '''
        with open(sample_sheet_path, 'r') as f:
            lines = f.readlines()
            lines = [x.strip() for x in lines]
            lines = [x.rstrip(',') for x in lines]
            lines = [x for x in lines if x]
            reads = []
            sentinel = False
            for i in range(0, len(lines)):
                if lines[i] == '[Reads]':
                    sentinel = True
                elif lines[i].startswith('['):
                    sentinel = False
                elif sentinel is True:
                    reads.append(lines[i])

            logging.debug("Reads: %s" % reads)

            if len(reads) == 1:
                return reads[0], None
            elif len(reads) == 2:
                return reads[0], reads[1]

            s = "_get_reads() found an unexpected number of values in {}"
            raise PipelineError(s.format(sample_sheet_path))
