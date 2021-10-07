from sequence_processing_pipeline.ConvertJob import ConvertJob
from sequence_processing_pipeline.QCJob import QCJob
from sequence_processing_pipeline.FastQCJob import FastQCJob
from sequence_processing_pipeline.GenPrepFileJob import GenPrepFileJob
from sequence_processing_pipeline.SequenceDirectory import SequenceDirectory
from sequence_processing_pipeline.PipelineError import PipelineError
import os
from os.path import join, exists
import json
from json.decoder import JSONDecodeError


class Pipeline:
    def __init__(self, configuration_file_path):
        try:
            f = open(configuration_file_path)
            self.configuration = json.load(f)['configuration']
            f.close()
        except FileNotFoundError:
            raise PipelineError(f'{configuration_file_path} does not exist.')
        except JSONDecodeError:
            raise PipelineError(f'{configuration_file_path} is not a valid ',
                                'json file')

        config = self.configuration['pipeline']
        self.younger_than = config['younger_than']
        self.older_than = config['older_than']
        self.archive_path = join(config['archive_path'])

        self._directory_check(self.archive_path, create=True)
        # self._directory_check(input_directory, create=False)
        # self._directory_check(output_directory, create=True)

        if self.older_than >= self.younger_than:
            raise PipelineError(
                'older_than cannot be equal to or less than younger_than.')

    def _directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    os.makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PipelineError, so it gets handled in the same location
                    # as the others.
                    raise PipelineError(str(e))
            else:
                raise PipelineError("directory_path '%s' does not exist." %
                                    directory_path)


if __name__ == '__main__':
    import logging

    configuration_file_path = ''
    # input_directory: A directory containing BCL files.
    input_directory = ''
    # A directory to store the Pipeline's products.
    output_directory = ''
    # identify Torque jobs using qiita_job_id
    qiita_job_id = ''
    # A path to the run_dir
    run_dir = ''
    # A path to the sample-sheet
    sample_sheet_path = ''
    # final_output_dir
    final_output_dir = ''

    # run_id='' might be useful for redefining run_id.
    # otherwise, it will be taken from the directory.
    sdo = SequenceDirectory(run_dir, sample_sheet_path, run_id='my_run_id')

    fastq_output_directory = sdo.fastq_results_directory

    pipeline = Pipeline(configuration_file_path)

    config = pipeline.configuration['bcl-convert']
    convert_job = ConvertJob(sdo.run_dir,
                             sdo.sample_sheet_path,
                             fastq_output_directory,
                             config['queue'],
                             config['nodes'],
                             config['nprocs'],
                             config['wallclock_time_in_hours'],
                             config['per_process_memory_limit'],
                             config['executable_path'],
                             config['modules_to_load'],
                             qiita_job_id)

    convert_job.run()

    config = pipeline.configuration['qc']
    qc_job = QCJob(sdo.run_dir,
                   sdo.sample_sheet_path,
                   config['mmi_db'],
                   final_output_dir,
                   config['queue'],
                   config['nodes'],
                   config['nprocs'],
                   config['wallclock_time_in_hours'],
                   config['per_process_memory_limit'],
                   config['fastp_executable_path'],
                   config['minimap2_executable_path'],
                   config['samtools_executable_path'],
                   config['modules_to_load'],
                   qiita_job_id)

    qc_job.run()

    # right now, it's preferable to use the trimmed/filtered products
    # as they are laid out on disk the way seqpro expects. Hence the
    # input directory for GenPrepFileJob is where QCJob moved its
    # output files to.
    #
    # seqpro needs a directory named RUN_ID w/subdirectories:
    # PROJECT_IDn/filtered_sequences or PROJECT_IDn/trimmed_sequences
    # to generate a prep-file for that project.
    gpf_input_path = join(final_output_dir, sdo.run_id)
    config = pipeline.configuration['seqpro']
    gpf_job = GenPrepFileJob(gpf_input_path,
                             sdo.sample_sheet_path,
                             final_output_dir,
                             config['seqpro_path'],
                             config['modules_to_load'],
                             qiita_job_id)

    gpf_job.run()

    config = pipeline.configuration['fastqc']
    output_directory = join(final_output_dir, 'FastQC')
    fastqc_job = FastQCJob(sdo.run_dir,
                           output_directory,
                           config['nprocs'],
                           config['nthreads'],
                           config['fastqc_executable_path'],
                           config['multiqc_executable_path'],
                           config['modules_to_load'],
                           qiita_job_id)

    fastqc_job.run()
