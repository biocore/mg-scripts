import logging
from os.path import exists, join
from sequence_processing_pipeline.QCCmdGenerator import QCCmdGenerator
from os import makedirs


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
            if self.h_filter is True:
                tmp = 'filtered_sequences'
            else:
                tmp = 'trimmed_sequences'

            if not exists(join(self.products_dir, tmp)):
                makedirs(join(self.products_dir, tmp), exist_ok=True)

        if self.a_trim is True:
            for fastq_file_path in self.fastq_file_paths:
                c = QCCmdGenerator(fastq_file_path, self.products_dir,
                                         self.project_name, self.nprocs,
                                         self.adapter_a,
                                         self.adapter_A, self.fastp_path,
                                         self.minimap2_path,
                                         self.samtools_path)

                if self.h_filter is True:
                    cmds.append(c.gen_chained_cmd(fastp_reports_dir,
                                                  self.human_phix_db_path))
                else:
                    cmds.append(c.gen_fastp_cmd())

            return cmds

        logging.warning("QCJob created w/a_trim set to False.")
        return None
