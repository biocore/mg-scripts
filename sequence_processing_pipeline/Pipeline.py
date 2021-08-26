import logging
from time import time as epoch_time
import os
import time
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory


class Pipeline:
    def __init__(self, scan_dir, threshold_in_hours=24):
        """
        Base class to define Pipelines for different Labs
        :param scan_dir:
        """
        logging.debug("Creating Pipeline Object")
        self.scan_dir = scan_dir
        self.rta_file_name = "RTAComplete.txt"
        logging.debug("Expected RTA File Name: %s" % self.rta_file_name)
        # internally, threshold will be represented in seconds
        self.threshold = threshold_in_hours * 60 * 60

    def _get_new_directories(self):
        """
        Scan for new sequencing raw data folders
        :return:
        """
        new_dirs = []
        current_time = epoch_time()

        for root, dirs, files in os.walk(self.scan_dir):
            for some_directory in dirs:
                some_path = os.path.join(root, some_directory)
                timestamp = os.path.getmtime(some_path)
                # if the last modified timestamp of the directory is within
                # the threshold of what is considered 'new data', then process
                # this directory.
                if((current_time - timestamp) < self.threshold):
                    formatted_ts = time.strftime('%m/%d/%Y :: %H:%M:%S', time.localtime(timestamp))
                    logging.debug("New directory found: %s\tTimestamp: %s" % (some_path, formatted_ts))
                    # save the path as well as the original epoch timestamp
                    # as a tuple.
                    if not os.path.exists(os.path.join(some_path, self.rta_file_name)):
                        logging.error("%s does not contain a file named '%s'." % (some_path, self.rta_file_name))
                    else:
                        logging.debug("Adding '%s' to processing list." % some_path)
                        new_dirs.append((some_path, timestamp))

        logging.debug("%d new directories found." % len(new_dirs))

        return new_dirs

    def process(self):
        """
        Process data.
        :return:
        """
        new_dirs = self._get_new_directories()

        for seq_dir, timestamp in new_dirs:
            sdo = SequenceDirectory(seq_dir, self.rta_file_name)
            if sdo.prepare_data_location():
                logging.debug("PREPARE DATA LOCATION RETURNED TRUE")

                results = sdo.process_data()
                '''

                # the idea behind this modification is to separate the act of
                # submitting the job from processing the data. Conceptually
                # they're two separate things. Also. it allows us an easy way to
                # test everything up to this point w/out submitting work to the
                # cluster. There may be other things we want to do with processed
                # data, for example.
                for result in results:
                    submission = SlurmBatch(result)
                    submission.submit()

                    # fastq_output = seq_dir
                    # human_filter(seq_dir, output_dir, fastq_output, nprocs)
                    # parse_csv?
                '''
            else:
                logging.debug("PREPARE DATA LOCATION RETURNED FALSE")
                '''
                # rely on print/log messages prior to this to explain error and notify users.
                # note that we can exit on the first sequence directory found that doesn't
                # pass prepare_data_location()'s sanity checks, or we could just note it and
                # keep going. Not sure what's the best approach.
                print("An error occured. Aborting...")
                exit(1)
                '''

