import logging
from os.path import exists, split, join  # , basename
from os import makedirs
from sequence_processing_pipeline.CmdGenerator import CmdGenerator


logging.basicConfig(level=logging.DEBUG)


#
# Be sure to have environment load modules before running the cmd string.
# module load fastp_0.20.1 samtools_1.12 minimap2_2.18' will be needed.
#


class fpmmp():
    def __init__(self, nprocs, trim_file, project_name, products_dir,
                 human_phix_db_path, adapter_a, adapter_A, a_trim, h_filter,
                 chemistry):
        '''
        Re-implements the non-copying functionality of fpmmp.py in Python.
        Generates an appropriate command-line string to run fastp and
        optionally minimap2 and samtools as needed.
        :param nprocs: Maximum number of processes/threads to use.
        :param trim_file: A file containing all of the R1 fastq.gz files.
        :param project_name: The name of the project. From sample-sheet.
        :param products_dir: The root directory to place products.
        :param human_phix_db_path: The path to human_phix_db_path.mmi
        :param adapter_a: Forward-read
        :param adapter_A: Reverse-read
        :param a_trim: A boolean. (Needs Adapter Trimming?)
        :param h_filter: A boolean. (Needs Human Filtering?)
        :param chemistry: Usually 'Default' or 'Amplicon'. From sample-sheet.
        '''

        # For now we'll assume that after using module, fastp, minimap2, and
        # samtools are all on the environment's path. We won't explicitly
        # give the full path to the commands in the result.

        self.project_name = project_name
        self.products_dir = products_dir
        self.adapter_a = adapter_a
        self.adapter_A = adapter_A
        self.a_trim = a_trim
        self.h_filter = h_filter
        self.chemistry = chemistry
        if nprocs > 16:
            nprocs = 16

        self.nprocs = str(nprocs)

        self.empty_file_list = None

        if adapter_a == 'NA':
            self.adapter_a = None

        if adapter_A == 'NA':
            self.adapter_A = None

        if not isinstance(self.a_trim, bool):
            raise ValueError("a_trim must be boolean.")

        if not isinstance(self.h_filter, bool):
            raise ValueError("h_filter must be boolean.")

        with open(trim_file, 'r') as f:
            lines = f.readlines()
            lines = [x.strip() for x in lines]
            self.trim_data = lines

        # added sanity check. Don't continue processing if the trim_file
        # contains a path to something unexpected.
        for line in self.trim_data:
            parent_dir, file_name = split(line)
            if file_name == '':
                # file_name will be '' if line ends in '/'
                raise ValueError('trim_file contains a directory.')

        # require path to human-phix-db.mmi even if it's not needed.
        # simpler to validate state.
        # possibly add checks to ensure file is valid
        self.human_phix_db_path = human_phix_db_path

        # there are two locations for empty_file_list. block2b() has 2 of its
        # four entries going to empty_file_list.txt in a different path:
        # CC: final_output renamed to products_dir
        # $final_output/filtered_sequences/${project}/empty_file_list.txt
        #
        # I think this is a mistake, and all four should go to below for
        # both block2a() and block2b().
        tmp = join(self.products_dir, self.project_name, 'empty_file_list.txt')
        self.empty_file_list_path = tmp

    def generate_commands(self):
        '''
        Generate the command-lines needed to QC the data, based on the
        parameters supplied to the object. The result will be the list of
        strings needed to process the fastq files referenced in the trim file.
        :return: A list of strings that can be run using Popen() and the like.
        '''
        logging.debug('QC.generate_command() called.')

        fastp_reports_dir = join(self.products_dir, 'fastp_reports_dir')

        # possible amplicon
        if self.a_trim is True:
            if self.h_filter is False:
                cmd = self._block2a(self.adapter_a,
                                    self.adapter_A,
                                    self.chemistry,
                                    self.project_name,
                                    self.products_dir)
            else:
                cmd = self._block2b(self.adapter_a,
                                    self.adapter_A,
                                    self.products_dir,
                                    self.project_name,
                                    fastp_reports_dir)
            return cmd

        logging.warning("QC.block2() called with a_trim == False.")
        return None

    def get_empty_file_list(self):
        '''
        Return the contents of the empty_file_list.txt file. Useful as an
        alternative to write_empty_file_list_to_path().
        :return: A list of empty files.
        '''
        return self.get_empty_file_list()

    def write_empty_file_list_to_path(self, path=None):
        '''
        Writes empty_file_list.txt to the supplied path.
        :param path: If path is None, the default path will be used.
        :return: the path to the written file
        '''
        some_path = path if path else self.empty_file_list_path

        if not self.empty_file_list:
            raise ValueError((f'There is nothing to write to {some_path}. '
                              'Run generate_command() first.'))

        with open(some_path, 'r') as f:
            for line in self.empty_file_list:
                f.write(f'{line}\n')

        return some_path

    def _block2a(self, adapter_a, adapter_A, chemistry, project_name,
                 products_dir):
        '''
        An internal method recreating the block of code in fpmmp.sh responsible
        for handling a_trim = True and h_filter = False.
        :param adapter_a: Forward-read
        :param adapter_A: Reverse-read
        :param chemistry: Chemistry of project. Usually 'Default' or 'Amplicon'
        :param project_name: Name of the project
        :param products_dir: The root directory to place products.
        :return: A list of strings to process the fastq files in trim_file.
        '''
        logging.debug('QC.block2a() called.')
        tmp = join(products_dir, project_name, 'trimmed_sequences')

        if not exists(tmp):
            logging.debug(f'creating {tmp}')
            makedirs(tmp, exist_ok=True)

        cmd_list = []
        # empty_file_list = []

        for fastq_file_path in self.trim_data:
            # technically, legacy will run fastp without adapter_a and
            # adapter_A if even one of then is None.
            # However, let's assume that both should be None and if
            # only one is None, then that's an Error. CmdGenerator will
            # test for that. Just pass the parameters.
            cmd_gen = CmdGenerator(fastq_file_path, products_dir,
                                   project_name, self.nprocs, adapter_a,
                                   adapter_A)
            cmd = cmd_gen.generate_fastp_cmd()
            logging.debug(f'Execute this command: {cmd}')
            cmd_list.append(cmd)

            '''
            partial = join(products_dir, project_name, 'trimmed_sequences')

            filename1_short, filename2_short = cmd.get_short_names()

            trimmed_file_path1 = filename1_short + '.fastp.fastq.gz'
            trimmed_file_path2 = filename2_short + '.fastp.fastq.gz'

            trimmed_sequence1 = join(partial, trimmed_file_path1)
            trimmed_sequence2 = join(partial, trimmed_file_path2)

            size1 = stat(trimmed_sequence1).st_size
            size2 = stat(trimmed_sequence2).st_size

            # amplicon/16s data. include index files
            if chemistry == 'Amplicon':
                # cp ${filename_index1}* $products_dir/${project}/
                #  trimmed_sequences
                # cp ${filename_index2}* $products_dir/${project}/
                #  trimmed_sequences
                pass
            else:
                logging.warning(
                    'QC.block2a() called w/chemistry != Amplicon')

            if size1 <= 500 or size2 <= 500:
                logging.warning(f'{filename1_short} is very small.')
                # mv $products_dir/${project}/trimmed_sequences/
                #  ${filename1_short}* $products_dir/${project}/
                #  trimmed_sequences/${filename2_short}*
                #   ${products_dir}/${project}/zero_files
                empty_file_list.append(trimmed_sequence1)
                empty_file_list.append(size1)
                empty_file_list.append(trimmed_sequence2)
                empty_file_list.append(size2)
            '''

        # self.empty_file_list = empty_file_list
        return cmd_list

    def _block2b(self, adapter_a, adapter_A, products_dir, project_name,
                 fastp_reports_dir):
        '''
        An internal method recreating the block of code in fpmmp.sh responsible
        for handling a_trim = True and h_filter = True.
        :param adapter_a: Forward-read
        :param adapter_A: Reverse-read
        :param products_dir: The root directory to place products.
        :param project_name: Name of the project
        :param fastp_reports_dir: The root directory to place fastp reports.
        :return: A list of strings to process the fastq files in trim_file.
        '''
        logging.debug('QC.block2b() called.')
        logging.debug("here i am inside true/true")

        tmp = join(products_dir, project_name, 'filtered_sequences')
        if not exists(tmp):
            makedirs(tmp, exist_ok=True)

        cmd_list = []
        empty_file_list = []

        for fastq_file_path in self.trim_data:
            cmd_gen = CmdGenerator(fastq_file_path, products_dir,
                                   project_name, self.nprocs, adapter_a,
                                   adapter_A)

            cmd = cmd_gen.generate_full_toolchain_cmd(fastp_reports_dir,
                                                      self.human_phix_db_path)

            # filename1 = basename(fastq_file_path)
            # filename1_short = basename(filename1) + '.fastq.gz'
            # filename2 = filename1.replace('_R1_00', '_R2_00')
            # filename2_short = basename(filename2) + '.fastq.gz'

            # cd parent_dir

            logging.debug(f'Execute this command: {cmd}')
            cmd_list.append(cmd)

            '''
            partial = join(products_dir, project_name, 'filtered_sequences')
            tmp1 = join(partial, filename1_short + '.trimmed.fastq.gz')
            size1 = stat(tmp1).st_size
            tmp2 = join(partial, filename2_short + '.trimmed.fastq.gz')
            size2 = stat(tmp2).st_size

            if size1 <= 500 or size2 <= 500:
                # mv something
                empty_file_list.append(tmp1)
                empty_file_list.append(size1)
                empty_file_list.append(tmp2)
                empty_file_list.append(size2)
            '''

        self.empty_file_list = empty_file_list

        return cmd_list


if __name__ == '__main__':
    logging.debug('test run started')

    nprocs = 16
    # assume current working directory is mg-scripts.
    trim_file = './sequence_processing_pipeline/tests/data/split_file_0'
    project_name = 'MyProject'
    products_dir = 'MyFinalOutputDir'
    human_phix_db_path = '/path/to/human-phix-db.mmi'
    adapter_a = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    adapter_A = 'GATCGGAAGAGCGTCGTGTAGGGAAAGGAGTGT'
    a_trim = True
    h_filter = True
    chemistry = 'Default'

    fpmmp = fpmmp(nprocs, trim_file, project_name, products_dir,
                  human_phix_db_path, adapter_a, adapter_A, a_trim, h_filter,
                  chemistry)

    cmds = fpmmp.generate_commands()
    for cmd in cmds:
        logging.debug(cmd)
    logging.debug('test run completed')
