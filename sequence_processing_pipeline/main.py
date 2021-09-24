import logging
from sequence_processing_pipeline.Pipeline import Pipeline


def main(input_directory, output_directory, final_output_directory,
         sample_sheet_path, younger_than, older_than, nprocs):
    '''
    Provide an easy entry point for users who would import the package and use
    it. It defines the parameters needed to process a directory of data.
    :param input_directory: A directory containing BCL files, or a directory
     of directories containing BCL files.
    :param output_directory: A directory to store the Pipeline's products.
    :param final_output_directory: A second directory to store products.
    :param younger_than: A threshold. Pipeline will not process directories
     younger than param hours.
    :param older_than: A threshold. Pipeline will not process directories
     older than param hours.
    :param nprocs: The number of parallel processes to run. Maximum value is
     currently 16.
    '''
    logging.debug('Sequence Processing Pipeline main() Started')

    pipeline = Pipeline(input_directory,
                        output_directory,
                        final_output_directory,
                        younger_than=younger_than,
                        older_than=older_than,
                        nprocs=nprocs)

    pipeline.process(sample_sheet_path)

    logging.debug('Sequence Processing Pipeline main() Completed')
