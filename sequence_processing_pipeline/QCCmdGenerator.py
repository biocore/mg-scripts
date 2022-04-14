from os.path import join, split
from os import makedirs


class QCCmdGenerator:
    def __init__(self, fastq_path, products_dir, project_name,
                 nprocs, adapter_a, adapter_A, fastp_path, minimap2_path,
                 samtools_path):
        """
        An object that generates the appropriate chain of commands, based on
        the parameters given to it.
        :param fastq_path: The path to a Fastq file.
        :param products_dir: The root directory to place products.
        :param project_name: The name of the project.
        :param nprocs: The maximum number of processes/threads to use.
        :param adapter_a: Forward-read
        :param adapter_A: Reverse-read
        """

        self.fastq_path = fastq_path
        self.current_dir = split(self.fastq_path)[0]
        self.parent_dir, self.filename1 = split(self.fastq_path)
        self.filename2 = self.filename1.replace('_R1_00', '_R2_00')
        self.products_dir = products_dir
        self.project_name = project_name
        self.nprocs = nprocs
        self.adapter_a = adapter_a
        self.adapter_A = adapter_A
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path

        if self.adapter_a is None:
            if self.adapter_A is not None:
                raise ValueError("adapter_a is declared but not adapter_A.")

        if self.adapter_A is None:
            if self.adapter_a is not None:
                raise ValueError("adapter_A is declared but not adapter_a.")

    def gen_fastp_cmd(self):
        """
        Generates a command-line string for running fastp, based on the
         parameters supplied to the object.
        :return: A string suitable for executing in Popen() and the like.
        """
        read1_input_path = join(self.current_dir, self.filename1)
        read2_input_path = join(self.current_dir, self.filename2)

        partial_path = join(self.current_dir, self.products_dir)
        json_output_path = join(partial_path, 'json',
                                self.filename1.replace('.fastq.gz', '.json'))
        html_output_path = join(partial_path, 'html',
                                self.filename1.replace('.fastq.gz', '.html'))
        read1_output_path = join(partial_path, 'trimmed_sequences',
                                 self.filename1.replace('.fastq.gz',
                                                        '.fastp.fastq.gz'))
        read2_output_path = join(partial_path, 'trimmed_sequences',
                                 self.filename2.replace('.fastq.gz',
                                                        '.fastp.fastq.gz'))
        report_title = self.filename1.replace('.fastq.gz', '') + '_report'

        result = self.fastp_path
        if self.adapter_a:
            # assume that if adapter_a is not None, then adapter_A is not None
            # as well. We performed this check already.
            result += (f' --adapter_sequence {self.adapter_a}'
                       f' --adapter_sequence_r2 {self.adapter_A}')

        result += (f' -l 100 -i {read1_input_path} -I {read2_input_path} -w '
                   f'{self.nprocs} -j {json_output_path} -h {html_output_path}'
                   f' -o {read1_output_path} -O {read2_output_path} -R '
                   f'{report_title}')

        return result

    def gen_chained_cmd(self, fastp_reports_dir,
                        human_phix_db_path):
        """
        Generates a command-line string for running fastp piping directly
         into minimap2 piping directly into samtools. The string is based
         on the parameters supplied to the object.
        :param fastp_reports_dir: Path to dir for storing json and html files
        :param human_phix_db_path: Path to the human_phix_db_path.mmi
        :return: A string suitable for executing in Popen() and the like.
        """

        read1_input_path = join(self.current_dir, self.filename1)
        read2_input_path = join(self.current_dir, self.filename2)

        tmp_path = join(self.project_name, fastp_reports_dir, 'json')
        makedirs(tmp_path, exist_ok=True)
        json_output_path = join(tmp_path, self.filename1.replace('.fastq.gz',
                                                                 '.json'))

        tmp_path = join(self.project_name, fastp_reports_dir, 'html')
        makedirs(tmp_path, exist_ok=True)
        html_output_path = join(tmp_path, self.filename1.replace('.fastq.gz',
                                                                 '.html'))

        partial = join(self.products_dir, 'filtered_sequences')

        path1 = join(partial, self.filename1.replace('.fastq.gz',
                                                     '.trimmed.fastq.gz'))
        path2 = join(partial, self.filename2.replace('.fastq.gz',
                                                     '.trimmed.fastq.gz'))

        result = self.fastp_path
        if self.adapter_a:
            # assume that if adapter_a is not None, then adapter_A is not None
            # as well. We performed this check already.
            result += (f' --adapter_sequence {self.adapter_a}'
                       f' --adapter_sequence_r2 {self.adapter_A}')

        result += (f' -l 100 -i {read1_input_path} -I {read2_input_path} '
                   f'-w {self.nprocs} -j {json_output_path} -h '
                   f'{html_output_path} --stdout | {self.minimap2_path} -ax'
                   f' sr -t {self.nprocs} {human_phix_db_path} - -a | '
                   f'{self.samtools_path} fastq -@ {self.nprocs} -f 12 -F '
                   f'256 -1 {path1} -2 {path2}')

        return result
