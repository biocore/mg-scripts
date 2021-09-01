import logging
from time import time as epoch_time
import os
import time
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from exceptions import PipelineError
from sequence_processing_pipeline.HumanFilter import HumanFilter


class Pipeline:
    def __init__(self, scan_dir, younger_than=48, older_than=24):
        """
        Base class to define Pipelines for different Labs
        :param scan_dir:
        """
        logging.debug("Creating Pipeline Object")
        self.scan_dir = scan_dir
        self.rta_file_name = "RTAComplete.txt"
        logging.debug("Expected RTA File Name: %s" % self.rta_file_name)
        # internally, threshold will be represented in seconds
        self.younger_than = younger_than * 60 * 60
        self.older_than = older_than * 60 * 60

    def _time_is_right(self, timestamp):
        # calculate how old the timestamp is
        delta_t = epoch_time() - timestamp

        # if the timestamp falls w/in the range defined by the object, then
        # it is a legitimate target.
        if delta_t < self.younger_than and delta_t > self.older_than:
            return True

        return False

    def _get_new_directories(self):
        """
        Scan for new sequencing raw data folders
        :return:
        """
        new_dirs = []

        for root, dirs, files in os.walk(self.scan_dir):
            for some_directory in dirs:
                some_path = os.path.join(root, some_directory)
                some_timestamp = os.path.getmtime(some_path)

                # whether a single BCL directory is passed, or a nested tree
                # of directories is passed, assume that a potential BCL
                # directory must contain a file named self.rta_file_name.
                if os.path.exists(os.path.join(some_path, self.rta_file_name)):
                    # if the last modified timestamp of the directory is
                    # within the threshold of what is considered 'new and
                    # completed data',then this directory is a legitimate
                    # target.
                    if self._time_is_right(some_timestamp):
                        formatted_ts = time.strftime('%m/%d/%Y %H:%M:%S (US/Pacific)', time.localtime(some_timestamp))
                        logging.info("Target found: %s\tTimestamp: %s" % (some_path, formatted_ts))
                        # save the path as well as the original epoch timestamp
                        # as a tuple.
                        new_dirs.append((some_path, some_timestamp))
                    else:
                        logging.debug("The timestamp for %s is not within bounds." % some_path)
                else:
                    # This is a warning, rather than an error because a
                    # directory of BCL directories would be a valid parameter,
                    # even though it doesn't contain BCL data itself.
                    logging.warning("%s does not contain a file named '%s'." % (some_path, self.rta_file_name))

        logging.debug("%d new directories found." % len(new_dirs))

        return new_dirs

    def process(self):
        """
        Process data.
        :return:
        """
        new_dirs = self._get_new_directories()

        for seq_dir, timestamp in new_dirs:
            try:
                sdo = SequenceDirectory(seq_dir, self.rta_file_name)
                sdo.prepare_data_location()
                sdo.process_data()
                hf = HumanFilter(sdo)
                hf.run()

                # this is technically it for the workflow itself.
                # the thing is, it's the human_filter and other? functions that do
                # sbatch and other Slurm calls to get the data going. We're going to
                # have some of this done by Qiita instead so it makes sense to just output
                # the slurm job array files and reorganize the code so that it makes the most
                # sense.

            except PipelineError as e:
                logging.error(e)
                # send out email notifications - or make the call to Qiita to let them know here.

