import logging
import shutil
from time import time as epoch_time
import os
import time
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.exceptions import PipelineError
from sequence_processing_pipeline.HumanFilter import HumanFilter
from sequence_processing_pipeline.util import system_call
from time import sleep
from datetime import datetime


class Pipeline:
    def __init__(self, root_dir, younger_than=48, older_than=24, should_filter=False):
        """
        Base class to define Pipelines for different Labs
        root_dir can be a single directory containing an RTAComplete.txt file
        or it can be a directory of directories containing RTAComplete.txt files.
        :param root_dir:
        """
        logging.debug("Creating Pipeline Object")
        self.root_dir = root_dir
        logging.debug("Root directory: %s" % self.root_dir)
        self.sentinel_file = "RTAComplete.txt"
        logging.debug("Sentinel File Name: %s" % self.sentinel_file)
        # internally, threshold will be represented in seconds
        self.younger_than = younger_than * 60 * 60
        self.older_than = older_than * 60 * 60
        logging.debug("Filter directories younger than %d hours old" % younger_than)
        logging.debug("Filter directories older than %d hours old" % older_than)
        self.should_filter = should_filter

    def _time_is_right(self, timestamp):
        # calculate how old the timestamp is
        delta_t = epoch_time() - timestamp

        # if the timestamp falls w/in the range defined by the object, then
        # it is a legitimate target.
        if delta_t < self.younger_than and delta_t > self.older_than:
            return True

        return False

    def _find_scan_directories(self):
        """
        Walk root directory and locate all scan directory root folders
        :return:
        """
        new_dirs = []

        for root, dirs, files in os.walk(self.root_dir):
            for some_directory in dirs:
                some_path = os.path.join(root, some_directory)
                if '/Data/Intensities/BaseCalls/' in some_path:
                    # the root directory of every scan directory will have
                    # this substring in the path of one or more of their
                    # subdirectories. By collecting these subdirectories
                    # and extracting the root directory of each one, we can
                    # build a set of unique root directories.
                    some_path = some_path.split('/Data/Intensities/BaseCalls')[0]
                    new_dirs.append(some_path)

        # remove duplicates
        new_dirs = list(set(new_dirs))

        for new_dir in new_dirs:
            logging.info("%s found." % new_dir)

        return new_dirs

    def _filter_directories_for_time(self, new_dirs):
        """
        Scan for new sequencing raw data folders
        :return:
        """
        filtered_dirs = []
        for some_path in new_dirs:
            # whether a single BCL directory is passed, or a nested tree
            # of directories is passed, assume that a potential BCL
            # directory must contain a file named self.sentinel_file.
            if os.path.exists(os.path.join(some_path, self.sentinel_file)):
                some_timestamp = os.path.getmtime(some_path)
                # if the last modified timestamp of the directory is
                # within the threshold of what is considered 'new and
                # completed data',then this directory is a legitimate
                # target.
                if self._time_is_right(some_timestamp):
                    formatted_ts = time.strftime('%m/%d/%Y %H:%M:%S (US/Pacific)', time.localtime(some_timestamp))
                    logging.info("Target found: %s\tTimestamp: %s" % (some_path, formatted_ts))
                    # save the path as well as the original epoch timestamp
                    # as a tuple.
                    filtered_dirs.append((some_path, some_timestamp))
                else:
                    logging.debug("The timestamp for %s is not within bounds." % some_path)
            else:
                # This is a warning, rather than an error because a
                # directory of BCL directories would be a valid parameter,
                # even though it doesn't contain BCL data itself.
                logging.warning("%s does not contain a file named '%s'." % (some_path, self.sentinel_file))

        logging.debug("%d new directories found." % len(filtered_dirs))

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
                elif sentinel == True:
                    # this should be a line in between
                    # [Bioinformatics] and the next section. Copy it
                    # to the buffer.
                    metadata.append(lines[i])
                # if the end of the file is reached before the next
                # header is found, that means there were no more
                # headers and that's okay too.

            # remove duplicate lines (this appears to be an issue in
            # the original bash scripts.)
            l = list(set(metadata))
            l.sort()

            # write the sanitized data out to legacy file.
            run_config_file_path = os.path.join(self.root_dir, 'run_config.txt')
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
                elif sentinel == True:
                    metadata.append(lines[i])

            l = list(set(metadata))
            l.sort()

            return metadata

    def get_date(self):
        '''
        Simulates the output of the Linux 'date' command.
        :return:
        '''
        # note the format of timestamp will be similar to:
        # 2021-09-01 21:17:19.005944
        # this is different than the format for 'date':
        # Wed Sep  1 21:11:17 PDT 2021
        # if this is an issue for consumers of this data, we can change it.
        timestamp = str(datetime.now())
        return timestamp

    def process(self):
        """
        Process data.
        :return:
        """
        new_dirs = self._find_scan_directories()
        new_dirs = self._filter_directories_for_time(new_dirs)

        for seq_dir, timestamp in new_dirs:
            try:
                sdo = SequenceDirectory(seq_dir)
                params = sdo.process_data()
                self._generate_run_config_file(params['sample_sheet_path'])
                contact_info = self._get_contact_information(params['sample_sheet_path'])
                self.prep_function()
                self.submit_job()

                #hf = HumanFilter(sdo)
                #hf.run()

                # this is technically it for the workflow itself.
                # the thing is, it's the human_filter and other? functions that do
                # sbatch and other Slurm calls to get the data going. We're going to
                # have some of this done by Qiita instead so it makes sense to just output
                # the slurm job array files and reorganize the code so that it makes the most
                # sense.

            except PipelineError as e:
                logging.error(e)
                # send out email notifications - or make the call to Qiita to let them know here.

    def prep_function(self):
        pass

    def _process_sacct_output(self, job_num, stdout):
        # treat stdout as a multiline string, and separate into a list of
        # strings.
        l = stdout.split('\n')
        # remove any leading or trailing whitespace in each line
        l = [x.strip() for x in l]
        # remove any line that doesn't begin with the job number.
        # there may be more than one row in the output:
        # 67767+0...
        # 67767+1...
        l = [x for x in l if x.startswith(job_num)]
        # convert each line into a list, so that we have a list of lists.
        # splitting on ' ' will create multiple empty strings due to multiple
        # spaces in between any two columns. Hence, use filter() to filter out
        # those empty strings using 'None' as a parameter.
        l = [list(filter(None, x.split(' '))) for x in l]

        return l

    def _was_job_successful(self, job_num, stdout):
        # process the stdout of sacct into something easy to work with
        results = self._process_sacct_output(job_num, stdout)
        # if one or more component jobs exited with a return code other than
        # 0, return False. Iterate through all jobs rather than return early
        # to write all jobs to logging.
        fail_flag = False
        for job in results:
            # the final element in each job list is the return code
            if job[-1] == '0:0':
                logging.debug("job %s exited with %s" % (str(job), job[-1]))
            else:
                logging.error("job %s exited with %s" % (str(job), job[-1]))
                fail_flag = True

        return fail_flag

    def submit_job(self, bcl_template, seq_dir, csvfile, base_mask, exp_name, bclconvert_template, polling_wait, root_sequence):
        if not os.path.exists(bcl_template):
            raise PipelineError("Can't locate job submission script %s (233)" % bcl_template)
        else:
            # assume this directory has already been created. we know it has been by this point.
            fastq_output = os.path.join(self.root_dir, 'Data', 'Fastq')
            if self.should_filter == False:
                cmd = ['sbatch',
                       '--parsable',
                       # presumably seq_proc is a defined quality of service level in the system
                       '--qos=seq_proc',
                       '--export=seqdir="%s",outputdir="%s",csvfile="%s",base_mask="%s"' % (
                       seq_dir, fastq_output, csvfile, base_mask),
                       '--job-name=%s' % exp_name,
                       '--partition=long %s' % bclconvert_template
                       ]

                stdout, stderr, return_code = system_call(cmd)
                # there may need some massaging to the output to turn it into a proper job_num
                job_num = stdout  # job_num should remain a string
                job_num_file = os.path.join(fastq_output, job_num) + '.txt'
                with open(job_num_file, 'w') as f:
                    f.write('start date == %s\n' % self.get_date())
                    f.write('experiment name == %s\n' % exp_name)
                    f.write('job_id == %s\n' % job_num)
                    f.write(
                        'submit args == sbatch --parsable --qos=seq_proc --export=seqdir="%s",outputdir="%s",csvfile="%s",base_mask="%s" --job-name=%s --partition=long %s' % (
                        seq_dir, fastq_output, csvfile, base_mask, exp_name, bclconvert_template))
                cmd = ['scontrol', 'show', 'job', job_num]
                while True:
                    stdout, stderr, return_code = system_call(cmd)
                    if return_code == -1:
                        break
                    else:
                        logging.debug(
                            "sbatch job %s still in progress. Sleeping %d seconds..." % (job_num, polling_wait))
                        sleep(polling_wait)
                with open(job_num_file, 'a') as f:
                    f.write('end date == %s\n' % self.get_date())

                cmd = ['saact', '-X', '-P', '-j', job_num]
                stdout, stderr, return_code = system_call(cmd)
                if self._was_job_successful(job_num, stdout):
                    logging.info("Initial processing complete for project %s %s" % (exp_name, root_sequence))
                    with open(job_num_file, 'a') as f:
                        # we could also do this for != 0 as well
                        f.write('exit status == 0 for %s\n' % job_num)
                    # this may freak out if the file is already there - we need to handle that.
                    # we should wrap this for known exceptions.
                    shutil.copy(csvfile, fastq_output)
                else:
                    logging.info("Initial processing issues %s %s %s" % (exp_name, root_sequence, csvfile))
        # TODO: Not sure yet if we should process processed and alockfile files here or
        #  at the end of process().
        #  touch ${seqdir}/processed && sleep 2
        #  rm ${seqdir}/alockfile