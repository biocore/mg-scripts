from sequence_processing_pipeline.BCLConvertJob import BCLConvertJob
from sequence_processing_pipeline.HumanFilterJob import HumanFilterJob
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.PipelineError import PipelineError
from time import time as epoch_time
import logging
import os
import pandas as pd


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
                job_params = self._prep_function(sample_sheet_params)

                bcl_convert_job = BCLConvertJob(self.root_dir,
                                                self.should_filter)

                bcl_convert_job.run(sample_sheet_params['sequence_directory'],
                                    sample_sheet_params['sample_sheet_path'],
                                    job_params['base_mask'],
                                    sample_sheet_params['experiment_name'],
                                    job_params['bclconvert_template'],
                                    job_params['root_sequence'])

                chemistry = sample_sheet_params['chemistry']
                slurm_array_task_id = 'UNKNOWN'

                human_filter_job = HumanFilterJob(sdo,
                                                  self.nprocs,
                                                  self.job_owner_home,
                                                  self.email_list,
                                                  chemistry,
                                                  self.output_dir,
                                                  self.final_output_dir)

                p1 = sample_sheet_params['sample_sheet_path']
                p2 = sample_sheet_params['sequence_directory']
                p3 = slurm_array_task_id
                human_filter_job.run(p1, p2, p3)

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
            job_params = self._prep_function(sample_sheet_params)

            bcl_convert_job = BCLConvertJob(self.root_dir,
                                            self.should_filter)

            bcl_convert_job.run(sample_sheet_params['sequence_directory'],
                                sample_sheet_params['sample_sheet_path'],
                                job_params['base_mask'],
                                sample_sheet_params['experiment_name'],
                                job_params['bclconvert_template'],
                                job_params['root_sequence'])

            chemistry = sample_sheet_params['chemistry']
            slurm_array_task_id = 'UNKNOWN'

            human_filter_job = HumanFilterJob(sdo,
                                              self.nprocs,
                                              self.job_owner_home,
                                              self.email_list,
                                              chemistry,
                                              self.output_dir,
                                              self.final_output_dir)

            human_filter_job.run(sample_sheet_path,
                                 sample_sheet_params['sequence_directory'],
                                 slurm_array_task_id)

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
                some_path = os.path.join(root, some_directory)
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
            if os.path.exists(os.path.join(some_path, self.sentinel_file)):
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
            s = os.path.join(self.root_dir, 'run_config.txt')
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

    def _process_sample_sheet(self, csv_file):
        '''

        :param csv_file:
        :return:
        '''
        logging.debug("Reading %s" % csv_file)
        df = pd.read_csv(csv_file, header=None)

        # drop all lines with no data (all "NaN"). Compresses sample sheet to
        # no blank lines.
        df = df.dropna(how='all')
        df = df.fillna(".")

        # Personally I feel initializing them to None is a better practice,
        # but we currently we rely on 0 to be the default value if the
        # elements are not found in the sheet.
        contact_index = 0
        bio_index = 0
        info_dict = {}

        # proposed default values for read_index 1 and 2, if values are not
        # found.
        read_index_1 = 151
        read_index_2 = 151

        for index in range(len(df)):
            if 'Bioinformatics' in df.iloc[index, 0]:
                bio_index = index
                print("bio_index=", index)
            elif '[Contact' in df.iloc[index, 0]:
                print("contact_index=", index)
                contact_index = index + 1
            else:
                contact_index = len(df)

            if '[Reads]' in df.iloc[index, 0]:
                read_index_1 = index
                read_index_2 = index + 3

            # assemble dictionary for csv variables

            keyword_list = ['Experiment', 'Assay', 'Chemistry',
                            'ReverseComplement']
            info_df = pd.DataFrame()

            for keyword in keyword_list:
                if keyword in df.iloc[index, 0]:
                    info_df = info_df.append(df.iloc[index])
                    value = df.iloc[index, 1]
                    info_dict[keyword] = value

        if info_dict["Chemistry"] == "Amplicon":
            logging.debug("Amplicon chemistry true. Removing false barcodes")

            n_count = None
            for line in csv_file:
                for entry in ["NNNNNNNN", "NNNNNNNNNNNN"]:
                    if entry in line:
                        n_count = str(len(entry))
                        logging.debug("N_count: %s" % n_count)
                        # assume once we find it in this file, we don't
                        # need to look further, as that would just
                        # overwrite the previous value.
                        break

            # just in case there are instances where n_count won't be found,
            # this preserves the original idea that 'I12' is the default.
            job_index_val = 'I12'
            map_n_count = {'12': 'I12', '8': 'I8', '0': 'I12'}

            if n_count:
                if n_count in map_n_count:
                    # if n_count found does not match a known value, we will
                    # default to 'I12'. If we should raise an Error instead,
                    # we should do that.
                    job_index_val = map_n_count[n_count]

            base_mask = "--use-bases-mask Y150,{},{},Y150"
            base_mask = base_mask.format(job_index_val, job_index_val)

            # check to see if both Read values are present at 150/151
            # if both, direction==2
            # else direction==1
            # used for bases-mask in bcl2fastq
            # pull from read_df or develop other method

            # replace N strings if Amplicon Sample sheets with blank for
            # processing.
            # can be moved to where we parse chemistry variable.
            # we have already copied original to backup location.
            # create new sample sheet with same name as original, removing
            # false barcodes
            csv_1 = open(csv_file, 'r')

            # should check to see if NNNNNNNN is present but this is historical
            csv_1 = ''.join([i for i in csv_1]).replace("NNNNNNNNNNNN", "")

            # I think instead of making a '.bak' file, reading the original
            # file for processing, and then reopening the file to be
            # overwritten, it would be better to simply read the original
            # file into memory (ala DF), close it, and write out the processed
            # contents of the original file to a new file. The names can
            # always be renamed afterward.
            csv_2 = open(csv_file, 'w')
            csv_2.writelines(csv_1)
            csv_2.close()

        bioinfo_df = pd.DataFrame(df.iloc[bio_index:contact_index - 1])

        logging.debug("contact_index == ", contact_index)

        contact_df = pd.DataFrame(df.iloc[contact_index:])
        if contact_df.empty:
            contact_df = pd.DataFrame({"some_email@gmail.com"})
        logging.debug(contact_df.empty)

        logging.debug("after empty df check")
        read_df = pd.DataFrame(df.iloc[read_index_1:read_index_2])

        logging.debug(bioinfo_df)
        logging.debug(contact_df)

        logging.debug(read_df.T)
        read_df = read_df.dropna(how='all')
        logging.debug(read_df)

    def _prep_function(self, sample_sheet_path):
        '''
        Get the metadata needed to submit a job.
        :param sample_sheet_path:
        :return:

        #rn1, rn2 = self._get_reads(sample_sheet_path)q
        #direction = 2 if rn2 else 1
        job_read_val = 'Y' + str(rn1)q
        column_count = self._get_column_count(sample_sheet_path)


          if [[ $column_count -eq "8" ]]; then
            if [[ $(awk '/NNNNNNNNNNNN/{getline;print;}' $csvfile | cut -f6 -d",") -eq "NNNNNNNNNNNN" ]]; then
              echo "strip line from file with false barcode"
            fi
            cutval="7"
          elif [[ $column_count > "8" ]]; then
            cutval=""
          fi

          index_value_6=$(awk '/Sample_ID/{getline;print;}' $csvfile | cut -f6 -d",")
          index_value_8=$(awk '/Sample_ID/{getline;print;}' $csvfile | cut -f8 -d",")
          index_size=${#index_value}
          index_6_size=${#index_value_6}

          if [[ -z $index_value_6 ]]; then
            job_index_val="I12"
          elif [[ $index_value_6 -eq "NNNNNNNNNNNN" ]]; then
            sed -i.bak s/NNNNNNNNNNNN//g $csvfile
            job_index_val="I12"
          elif [[ $index_6_size -eq "8" ]]; then
            job_index_val="I8"
          elif [[ $index_6_size -eq "12" ]]; then
            job_index_val="I12"
          fi

          # if index 7 == NNNNNNNNNNNN then must be metagenomic with multiple rows
          # because of additional first column
          if [[ $index_value_7 -eq "NNNNNNNNNNNN" ]]; then
            sed -i.bak s/NNNNNNNNNNNN/g $csvfile
            job_index_val="I12"
          fi

          if [[ $direction -eq "1" ]]; then
            base_mask="--use-bases-mask Y$rn1,$job_index_val"
          elif [[ $direction -eq "2" ]]; then
            base_mask="--use-bases-mask Y$rn1,$job_index_val,$job_index_val,Y$rn1"
          fi

          chemistry=$(awk -F',' '/Chemistry/{print $2}' $csvfile)
          echo chemistry==$chemistry

          root_sequence=$(basename $dirpath)

          email_list=null
          job_o_out=localhost:/home/jede9131/seq_jobs/$(basename ${seqdir})
          job_e_out=localhost:/home/jede9131/seq_jobs/$(basename ${seqdir})



        # revisit.
        index_value_6 = "NOTHING"
        index_value_7 = "NOTHING"
        root_sequence = 'NOTHING'

        if index_value_6 == None:
            job_index_val = 'I12'
        elif index_value_6 == "NNNNNNNNNNNN":
            # remove string from csvfile
            # sed -i.bak s/NNNNNNNNNNNN//g $csvfile
            job_index_val = "I12"
        elif index_value_6 == 8:
            job_index_val = "I8"
        elif index_value_6 == 12:
            job_index_val = "I12"

        if index_value_7 == "NNNNNNNNNNNN":
            # sed -i.bak s/NNNNNNNNNNNN/g $csvfile
            job_index_val = "I12"

        # assuming direction can only equal 1 or 2:
        base_mask = "--use-bases-mask Y%s,%s" % (rn1, job_index_val) if direction == 1 else "--use-bases-mask Y%s,%s,%s,Y%s" % (rn1, job_index_val, job_index_val, rn1)

        # this script may need to be obtained from the admin
        bclconvert_template = "/home/somebody/seq_proc_dev/bclconvert_human_slurm.sh"

        job_params = { 'base_mask': base_mask, 'bclconvert_template': bclconvert_template, 'root_sequence': root_sequence }

        return job_params
        '''
        pass

