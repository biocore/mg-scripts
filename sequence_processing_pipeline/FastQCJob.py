from multiprocessing import Pool
import logging
from os import makedirs, listdir
from os.path import join
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.exec_fastqc import exec_fastqc
from sequence_processing_pipeline.exec_multiqc import exec_multiqc


class FastQCJob(Job):
    def __init__(self, run_dir, output_directory, nprocs, nthreads,
                 fastqc_path, multiqc_path, modules_to_load):
        self.job_name = 'FastQCJob'
        super().__init__(run_dir, self.job_name, [fastqc_path, multiqc_path],
                         modules_to_load)
        self.nprocs = str(nprocs)
        self.nthreads = str(nthreads)
        self.fastqc_path = fastqc_path
        self.multiqc_path = multiqc_path
        self.modules_to_load = modules_to_load
        # assume that the path to raw fastqs is run_dir/Data/Fastq
        self.raw_fastq_path = join(self.run_dir, 'Data', 'Fastq')
        # output_directory is the root path for creating a FastQC directory
        # with FastQC output. It's also the location for finding trimmed and
        # filtered fastq files as well.
        self.fastqc_output_path = output_directory
        # '.../output_dir/run_id/project/filter_type'
        self.processed_fastq_path = output_directory

    def run(self):

        # gather the parameters for processing all relevant raw fastq files.
        params = self._scan_fastq_files(self.raw_fastq_path,
                                        self.fastqc_output_path,
                                        is_raw_input=True)

        # next, do the same for the trimmed/filtered fastq files.
        params += self._scan_fastq_files(self.processed_fastq_path,
                                         self.fastqc_output_path)

        # Pool emulates the behavior of the GNU parallel command that was
        # used in the legacy scripts.
        with Pool(processes=self.nprocs) as fastqc_pool:
            fastqc_pool.map_async(exec_fastqc, params).get()

        # For now, create an entirely separate pool for multiqc from fastqc
        # to ensure all fastqc calls are completed before multiqc ones start.

        # gather the parameters for processing all relevant raw fastq files.
        params = self._scan_fastq_files(self.raw_fastq_path,
                                        self.fastqc_output_path,
                                        is_raw_input=True,
                                        for_multiqc=True)

        # next, do the same for the trimmed/filtered fastq files.
        params += self._scan_fastq_files(self.processed_fastq_path,
                                         self.fastqc_output_path,
                                         for_multiqc=True)

        with Pool(processes=self.nprocs) as multiqc_pool:
            multiqc_pool.map_async(exec_multiqc, params).get()

    def _find_projects(self, path_to_run_id_data_fastq_dir):
        results = []
        for directory in listdir(path_to_run_id_data_fastq_dir):
            project_dir = join(path_to_run_id_data_fastq_dir, directory)
            files = self._find_files(project_dir)

            # extract only fastq files from the list
            files = [x for x in files if x.endswith('.fastq.gz')]

            if files:
                tmp = ' '.join(files)
                if 'trimmed_sequences' in tmp:
                    # results are a_trim = True, h_filter= = False
                    filter_type = 'trimmed_sequences'
                elif 'filtered_sequences' in tmp:
                    # results are a_trim = True, h_filter= = True
                    # (trimmed AND filtered)
                    filter_type = 'filtered_sequences'
                elif 'amplicon' in tmp:
                    # results are a_trim = False, h_filter= = False
                    filter_type = 'amplicon'
                else:
                    raise ValueError("indeterminate type")

                if filter_type != 'amplicon':
                    # filter out index '_In_' files
                    files = [x for x in files if '_R1_' in x or '_R2_' in x]

                logging.debug("%s is a project directory" % project_dir)
                results.append((directory, filter_type, project_dir, files))

        return results

    def _scan_fastq_files(self, input_path, fastqc_output_path,
                          is_raw_input=False, for_multiqc=False):
        '''
        Scan a .../run_id/Data/Fastq directory or an .../output_dir/run_id
        directory for FastQC job information.
        :param input_path: The path to a directory w/Fastq files.
        :param fastqc_output_path: The path to store all FastQC output.
        :param is_raw_input: Boolean to modify the output path.
        :param for_multiqc: Boolean. Generate params for FastQC or MultiQC.
        :return: A list of tuples suitable for adding to a subprocess.Pool
        '''
        projects = self._find_projects(input_path)

        fastqc_results = []
        multiqc_results = []

        for project_name, filter_type, fastq_path, files in projects:
            tmp = join(fastqc_output_path, project_name)
            if is_raw_input:
                output_dir = join(tmp, 'bclconvert')
            else:
                output_dir = join(tmp, filter_type)
            makedirs(output_dir, exist_ok=True)
            for some_file in files:
                if for_multiqc:
                    fastq_trimmed_dir = join(self.processed_fastq_path,
                                             project_name,
                                             filter_type)

                    multiqc_results.append((self.nthreads,
                                            some_file,
                                            fastqc_output_path,
                                            project_name,
                                            filter_type,
                                            fastq_trimmed_dir))
                else:
                    fastqc_results.append((self.nthreads, some_file,
                                           output_dir))

        return fastqc_results, multiqc_results
