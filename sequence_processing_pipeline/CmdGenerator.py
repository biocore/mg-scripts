from os.path import join, split


class CmdGenerator:
    def __init__(self, fastq_path, products_dir, project_name,
                 nprocs, adapter_a, adapter_A):
        '''

        :param fastq_path: The path to a Fastq file.
        :param products_dir: The root directory to place products.
        :param project_name: The name of the project.
        :param nprocs: The maximum number of processes/threads to use.
        :param adapter_a: Forward-read
        :param adapter_A: Reverse-read
        '''

        self.fastq_path = fastq_path
        self.current_dir = split(self.fastq_path)[0]
        self.parent_dir, self.filename1 = split(self.fastq_path)
        self.filename2 = self.filename1.replace('_R1_00', '_R2_00')
        self.products_dir = products_dir
        self.project_name = project_name
        self.nprocs = nprocs
        self.adapter_a = adapter_a
        self.adapter_A = adapter_A
        self.filename1_short = self.filename1.replace('.fastq.gz', '')
        self.filename2_short = self.filename2.replace('.fastq.gz', '')

        if self.adapter_a is None:
            if self.adapter_A is not None:
                raise ValueError("adapter_a is declared but not adapter_A.")

        if self.adapter_A is None:
            if self.adapter_a is not None:
                raise ValueError("adapter_A is declared but not adapter_a.")

    def get_short_names(self):
        '''
        Return the shortnames for the given Fastq file. Used in a number of
         places so it's preferable to obtain them from a single source.
        :return: a tuple of (filename1_short, filename2_short)
        '''
        return (self.filename1_short, self.filename2_short)

    def generate_fastp_cmd(self):
        '''
        Generates a command-line string for running fastp, based on the
         parameters supplied to the object.
        :return: A string suitable for executing in Popen() and the like.
        '''
        read1_input_path = join(self.current_dir, self.filename1)
        read2_input_path = join(self.current_dir, self.filename2)

        partial_path = join(self.current_dir, self.products_dir,
                            self.project_name)
        json_output_path = join(partial_path, 'json',
                                self.filename1_short + '.json')
        html_output_path = join(partial_path, 'html',
                                self.filename1_short + '.html')
        read1_output_path = join(partial_path, 'trimmed_sequences',
                                 self.filename1_short + '.fastp.fastq.gz')
        read2_output_path = join(partial_path, 'trimmed_sequences',
                                 self.filename2_short + '.fastp.fastq.gz')
        report_title = self.filename1_short + '_report'

        self.fastp_cmd_list = ['fastp']
        if self.adapter_a:
            # assume that if adapter_a is not None, then adapter_A is not None
            # as well. We performed this check already.
            self.fastp_cmd_list.append('--adapter_sequence')
            self.fastp_cmd_list.append(self.adapter_a)
            self.fastp_cmd_list.append('--adapter_sequence_r2')
            self.fastp_cmd_list.append(self.adapter_A)

        self.fastp_cmd_list.append('-l')
        self.fastp_cmd_list.append('100')
        self.fastp_cmd_list.append('-i')
        self.fastp_cmd_list.append(read1_input_path)
        self.fastp_cmd_list.append('-I')
        self.fastp_cmd_list.append(read2_input_path)
        self.fastp_cmd_list.append('-w')
        self.fastp_cmd_list.append(self.nprocs)
        self.fastp_cmd_list.append('-j')
        self.fastp_cmd_list.append(json_output_path)
        self.fastp_cmd_list.append('-h')
        self.fastp_cmd_list.append(html_output_path)
        self.fastp_cmd_list.append('-o')
        self.fastp_cmd_list.append(read1_output_path)
        self.fastp_cmd_list.append('-O')
        self.fastp_cmd_list.append(read2_output_path)
        self.fastp_cmd_list.append('-R')
        self.fastp_cmd_list.append(report_title)

        return ' '.join(self.fastp_cmd_list)

    def generate_full_toolchain_cmd(self, fastp_reports_dir,
                                    human_phix_db_path):
        '''
        Generates a command-line string for running fastp piping directly
         into minimap2 piping directly into samtools. The string is based
         on the parameters supplied to the object.
        :param fastp_reports_dir: Path to dir for storing json and html files
        :param human_phix_db_path: Path to the human_phix_db_path.mmi
        :return: A string suitable for executing in Popen() and the like.
        '''

        read1_input_path = join(self.current_dir, self.filename1)
        read2_input_path = join(self.current_dir, self.filename2)

        json_output_path = join(fastp_reports_dir, self.project_name, 'json',
                                self.filename1_short + '.json')
        html_output_path = join(fastp_reports_dir, self.project_name, 'html',
                                self.filename1_short + '.html')

        partial = join(self.products_dir, self.project_name,
                       'filtered_sequences')

        # TODO: this string might be different. confirm it needs
        #  trimmed.fastq.gz
        path1 = join(partial, self.filename1_short + '.trimmed.fastq.gz')
        path2 = join(partial, self.filename2_short + '.trimmed.fastq.gz')

        self.fastp_cmd_list = ['fastp']
        if self.adapter_a:
            # assume that if adapter_a is not None, then adapter_A is not None
            # as well. We performed this check already.
            self.fastp_cmd_list.append('--adapter_sequence')
            self.fastp_cmd_list.append(self.adapter_a)
            self.fastp_cmd_list.append('--adapter_sequence_r2')
            self.fastp_cmd_list.append(self.adapter_A)
        self.fastp_cmd_list.append('-l')
        self.fastp_cmd_list.append('100')
        self.fastp_cmd_list.append('-i')
        self.fastp_cmd_list.append(read1_input_path)
        self.fastp_cmd_list.append('-I')
        self.fastp_cmd_list.append(read2_input_path)
        self.fastp_cmd_list.append('-w')
        self.fastp_cmd_list.append(self.nprocs)
        self.fastp_cmd_list.append('-j')
        self.fastp_cmd_list.append(json_output_path)
        self.fastp_cmd_list.append('-h')
        self.fastp_cmd_list.append(html_output_path)
        self.fastp_cmd_list.append('--stdout')

        self.minimap_cmd_list = ['minimap2']
        self.minimap_cmd_list.append('-ax')
        self.minimap_cmd_list.append('sr')
        self.minimap_cmd_list.append('-t')
        self.minimap_cmd_list.append(self.nprocs)
        self.minimap_cmd_list.append(human_phix_db_path)
        self.minimap_cmd_list.append('-')
        self.minimap_cmd_list.append('-a')

        self.samtools_cmd_list = ['samtools']
        self.samtools_cmd_list.append('fastq')
        self.samtools_cmd_list.append('-@')
        self.samtools_cmd_list.append(self.nprocs)
        self.samtools_cmd_list.append('-f')
        self.samtools_cmd_list.append('12')
        self.samtools_cmd_list.append('-F')
        self.samtools_cmd_list.append('256')
        self.samtools_cmd_list.append('-1')
        self.samtools_cmd_list.append(path1)
        self.samtools_cmd_list.append('-2')
        self.samtools_cmd_list.append(path2)

        # create the final command piping all other commands together.
        return ' '.join(self.fastp_cmd_list) + ' | ' + ' '.join(
            self.minimap_cmd_list) + ' | ' + ' '.join(self.samtools_cmd_list)
