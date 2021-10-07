from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os.path import basename
from sequence_processing_pipeline.PipelineError import PipelineError
from time import time as epoch_time
import logging
import os
from os.path import join, exists


class SequenceDirectory:
    def __init__(self, run_dir, sample_sheet_path, run_id=None):
        self.sentinel_file = "RTAComplete.txt"

        if run_dir:
            if exists(run_dir):
                self.run_dir = run_dir
            else:
                raise PipelineError("Directory %s does not exist." % run_dir)
        else:
            raise PipelineError("A value for run_dir must be supplied.")

        if run_id:
            # potentially useful for run_dirs that aren't following naming
            # convention.
            self.run_id = run_id
        else:
            self.run_id = basename(self.run_dir)

        if sample_sheet_path:
            if exists(sample_sheet_path):
                self.sample_sheet_path = sample_sheet_path
            else:
                raise PipelineError(f"External sample sheet "
                                    f"{sample_sheet_path} does not exist.")
        else:
            raise PipelineError("An external sample sheet must be supplied.")

        s = join(self.run_dir, 'Data', 'Fastq')
        self.fastq_results_directory = s

        try:
            os.makedirs(self.fastq_results_directory, exist_ok=True)
        except OSError as e:
            # this is a known potential error. Re-raise it as a
            # PipelineError, so it gets handled in the same location as the
            # others.
            raise PipelineError(str(e))

        # extract needed metadata from the sample sheet.
        sheet = KLSampleSheet(self.sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            raise PipelineError(f"Sample sheet {self.sample_sheet_path} "
                                "is not valid.")

        header = valid_sheet.Header
        self.experiment_name = header['Experiment Name']
        self.chemistry = header['Chemistry']

    def find_bcl_directories(self):
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

        return new_dirs

    def _filter_directories_for_time(self, new_dirs, younger_than, older_than):
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
                delta_t = epoch_time() - some_timestamp
                if delta_t < younger_than and delta_t > older_than:
                    # save the path as well as the original epoch timestamp
                    # as a tuple.
                    filtered_dirs.append((some_path, some_timestamp))
                else:
                    logging.debug(f"The timestamp for {some_path} is not"
                                  "within bounds.")
            else:
                # This is a warning, rather than an error because a
                # directory of BCL directories would be a valid parameter,
                # even though it doesn't contain BCL data itself.
                s = "{} does not contain a file named '{}'."
                logging.warning(s.format(some_path, self.sentinel_file))

        return filtered_dirs
