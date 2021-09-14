import logging
import sys
from sequence_processing_pipeline.KnightLabPipeline import KnightLabPipeline


def main(input_directory, output_directory, final_output_directory, younger_than, older_than, nprocs):
    '''
    Provide an easy entry point for users who would import the package and use
    it. It defines the parameters needed to process a directory of data.
    :param input_directory: A directory containing BCL files, or a directory
     of directories containing BCL files.
    :param output_directory: A directory to store the Pipeline's products.
    :param final_output_directory: A second directory to store the Pipeline's products.
    :param younger_than: A threshold. Pipeline will not process directories
     younger than param hours.
    :param older_than: A threshold. Pipeline will not process directories
     older than param hours.
    :param nprocs: The number of parallel processes to run. Maximum value is
     currently 16.
    '''

    logging.debug('Sequence Processing Pipeline Main() Started')

    pipeline = KnightLabPipeline(input_directory,
                                 output_directory,
                                 final_output_directory,
                                 younger_than=younger_than,
                                 older_than=older_than,
                                 nprocs=nprocs)
    pipeline.process()

    logging.debug('Sequence Processing Pipeline Main() Completed')


def single_directory_main(input_directory, output_directory, final_output_directory,
                          younger_than, older_than, nprocs):
    '''
    Provide an easy entry point for users who would import the package and use
    it. It defines the parameters needed to process a directory of data.
    :param input_directory: A directory containing BCL files, or a directory
     of directories containing BCL files.
    :param output_directory: A directory to store the Pipeline's products.
    :param final_output_directory: A second directory to store the Pipeline's products.
    :param younger_than: A threshold. Pipeline will not process directories
     younger than param hours.
    :param older_than: A threshold. Pipeline will not process directories
     older than param hours.
    :param nprocs: The number of parallel processes to run. Maximum value is
     currently 16.
    '''
    logging.debug('Sequence Processing Pipeline SingleDirectoryMain() Started')

    pipeline = KnightLabPipeline(input_directory,
                                 output_directory,
                                 final_output_directory,
                                 younger_than=younger_than,
                                 older_than=older_than,
                                 nprocs=nprocs)
    pipeline.process_one()

    logging.debug('Sequence Processing Pipeline SingleDirectoryMain() Completed')

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    # logging.basicConfig(level=logging.DEBUG, filename='sequence_processing_pipeline.log')

    main(sys.argv[1], sys.argv[2], sys.argv[3], 48, 24, 16)
