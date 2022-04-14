import logging
from os.path import exists
from os import makedirs
from os.path import join, split


class QCHelper():
    def __init__(self, nprocs, fastq_file_paths, project_name, products_dir,
                 human_phix_db_path, adapter_a, adapter_A, a_trim, h_filter,
                 chemistry, fastp_path, minimap2_path, samtools_path):
        """
        Re-implements the non-copying functionality of fpmmp.py in Python.
        Generates an appropriate command-line string to run fastp and
        optionally minimap2 and samtools as needed.
        :param nprocs: Maximum number of processes/threads to use.
        :param fastq_file_paths: A list of paths to the R1 fastq.gz files.
        :param project_name: The name of the project. From sample-sheet.
        :param products_dir: The root directory to place products.
        :param human_phix_db_path: The path to human_phix_db_path.mmi
        :param adapter_a: Forward-read
        :param adapter_A: Reverse-read
        :param a_trim: A boolean. (Needs Adapter Trimming?)
        :param h_filter: A boolean. (Needs Human Filtering?)
        :param chemistry: Usually 'Default' or 'Amplicon'. From sample-sheet.
        :param fastp_path: The path to the fastp executable
        :param minimap2_path: The path to the minimap2 executable
        :param samtools_path: The path to the samtools executable
        """

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
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path

        if nprocs > 16:
            nprocs = 16

        self.nprocs = str(nprocs)

        if adapter_a == 'NA':
            self.adapter_a = None

        if adapter_A == 'NA':
            self.adapter_A = None

        if self.adapter_a is None:
            if self.adapter_A is not None:
                raise ValueError("adapter_a is declared but not adapter_A.")

        if self.adapter_A is None:
            if self.adapter_a is not None:
                raise ValueError("adapter_A is declared but not adapter_a.")

        if not isinstance(self.a_trim, bool):
            raise ValueError("a_trim must be boolean.")

        if not isinstance(self.h_filter, bool):
            raise ValueError("h_filter must be boolean.")

        self.fastq_file_paths = [x.strip() for x in fastq_file_paths if
                                 '_R1_' in x]

        # require path to human-phix-db.mmi even if it's not needed.
        # simpler to validate state.
        # possibly add checks to ensure file is valid
        self.human_phix_db_path = human_phix_db_path

        if self.h_filter is True:
            self.outdir_name = 'filtered_sequences'
        else:
            self.outdir_name = 'trimmed_sequences'

    def generate_commands(self):
        """
        Generate the command-lines needed to QC the data, based on the
        parameters supplied to the object. The result will be the list of
        strings needed to process the fastq files referenced in the trim file.
        :return: A list of strings that can be run using Popen() and the like.
        """
        fastp_reports_dir = join(self.products_dir, 'fastp_reports_dir')

        # if this directory is not made, then fastp will not create the html
        # and json directories and generate output files for them.
        makedirs(fastp_reports_dir, exist_ok=True)

        cmds = []

        if self.a_trim is True:
            if not exists(join(self.products_dir, self.outdir_name)):
                makedirs(join(self.products_dir, self.outdir_name),
                         exist_ok=True)

        if self.a_trim is True:
            for fastq_file_path in self.fastq_file_paths:
                current_dir = split(fastq_file_path)[0]
                _, filename1 = split(fastq_file_path)
                filename2 = filename1.replace('_R1_00', '_R2_00')

                if self.h_filter is True:
                    cmds.append(self._gen_chained_cmd(current_dir,
                                                      filename1,
                                                      filename2,
                                                      fastp_reports_dir))
                else:
                    cmds.append(self._gen_fastp_cmd(current_dir,
                                                    filename1,
                                                    filename2))

            return cmds

        logging.warning("QCJob created w/a_trim set to False.")
        return None

    def _gen_fastp_cmd(self, current_dir, filename1, filename2):
        """
        Generates a command-line string for running fastp, based on the
         parameters supplied to the object.
        :return: A string suitable for executing in Popen() and the like.
        """
        read1_input_path = join(current_dir, filename1)
        read2_input_path = join(current_dir, filename2)

        partial_path = join(current_dir, self.products_dir)
        json_output_path = join(partial_path, 'json',
                                filename1.replace('.fastq.gz', '.json'))
        html_output_path = join(partial_path, 'html',
                                filename1.replace('.fastq.gz', '.html'))
        read1_output_path = join(partial_path, 'trimmed_sequences',
                                 filename1.replace('.fastq.gz',
                                                   '.fastp.fastq.gz'))
        read2_output_path = join(partial_path, 'trimmed_sequences',
                                 filename2.replace('.fastq.gz',
                                                   '.fastp.fastq.gz'))
        report_title = filename1.replace('.fastq.gz', '') + '_report'

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

    def _gen_chained_cmd(self, current_dir, filename1, filename2,
                         fastp_reports_dir):
        """
        Generates a command-line string for running fastp piping directly
         into minimap2 piping directly into samtools. The string is based
         on the parameters supplied to the object.
        :param fastp_reports_dir: Path to dir for storing json and html files
        :param human_phix_db_path: Path to the human_phix_db_path.mmi
        :return: A string suitable for executing in Popen() and the like.
        """

        read1_input_path = join(current_dir, filename1)
        read2_input_path = join(current_dir, filename2)

        tmp_path = join(self.project_name, fastp_reports_dir, 'json')
        makedirs(tmp_path, exist_ok=True)
        json_output_path = join(tmp_path, filename1.replace('.fastq.gz',
                                                            '.json'))

        tmp_path = join(self.project_name, fastp_reports_dir, 'html')
        makedirs(tmp_path, exist_ok=True)
        html_output_path = join(tmp_path, filename1.replace('.fastq.gz',
                                                            '.html'))

        partial = join(self.products_dir, 'filtered_sequences')

        path1 = join(partial, filename1.replace('.fastq.gz',
                                                '.trimmed.fastq.gz'))
        path2 = join(partial, filename2.replace('.fastq.gz',
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
                   f' sr -t {self.nprocs} {self.human_phix_db_path} - -a | '
                   f'{self.samtools_path} fastq -@ {self.nprocs} -f 12 -F '
                   f'256 -1 {path1} -2 {path2}')

        return result
