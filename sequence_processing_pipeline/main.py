import logging
import sys
from sequence_processing_pipeline.KnightLabPipeline import KnightLabPipeline


def main(input_directory, output_directory, younger_than, older_than, nprocs,
         should_filter):
    '''
    Provide an easy entry point for users who would import the package and use
    it. It defines the parameters needed to process a directory of data.
    :param input_directory: A directory containing BCL files, or a directory
     of directories containing BCL files.
    :param output_directory: A directory to store the Pipeline's products.
    :param younger_than: A threshold. Pipeline will not process directories
     younger than param hours.
    :param older_than: A threshold. Pipeline will not process directories
     older than param hours.
    :param nprocs: The number of parallel processes to run. Maximum value is
     currently 16.
    :param should_filter: TODO: Not sure why this shouldn't always be True.
    :return: TODO: Could return True/False if the processing was successful.
      However, logging is more substantial. Exceptions are better for fatal
      errors. Hard to say if anything should be returned. If it's fatal it
      should be returned as an exception. If it's not fatal, then it should
      be logged because there can be numerous non-fatal warnings that
      happen during a run.
    '''

    # TODO: We should have a parameter should_contain and have it be a list of
    #  strings that represent filenames used as sentinels.
    # TODO: How do we handle passing sample sheets? We need to pass them
    #  through this main, but we also need to handle the case where we have
    #  multiple directories. Perhaps if we need to rearrange directory of
    #  directories so that it's above this main. and then we add a param to
    #  supply the path to the sample file, which can be inside the directory
    #  or external. We could also have a switch to say, file is in the
    #  directory.
    logging.debug('Sequence Processing Pipeline Main() Started')

    pipeline = KnightLabPipeline(input_directory,
                                 output_directory,
                                 younger_than=younger_than,
                                 older_than=older_than,
                                 should_filter=should_filter,
                                 nprocs=nprocs)
    pipeline.process()

    logging.debug('Sequence Processing Pipeline Main() Completed')


def single_directory_main(input_directory, output_directory, younger_than,
                          older_than, nprocs, should_filter):
    '''
    Provide an easy entry point for users who would import the package and use
    it. It defines the parameters needed to process a directory of data.
    :param input_directory: A directory containing BCL files, or a directory
     of directories containing BCL files.
    :param output_directory: A directory to store the Pipeline's products.
    :param younger_than: A threshold. Pipeline will not process directories
     younger than param hours.
    :param older_than: A threshold. Pipeline will not process directories
     older than param hours.
    :param nprocs: The number of parallel processes to run. Maximum value is
     currently 16.
    :param should_filter: TODO: Not sure why this shouldn't always be True.
    :return: TODO: Could return True/False if the processing was successful.
      However, logging is more substantial. Exceptions are better for fatal
      errors. Hard to say if anything should be returned. If it's fatal it
      should be returned as an exception. If it's not fatal, then it should
      be logged because there can be numerous non-fatal warnings that
      happen during a run.
    '''

    # TODO: We should have a parameter should_contain and have it be a list of
    #  strings that represent filenames used as sentinels.
    # TODO: How do we handle passing sample sheets? We need to pass them
    #  through this main, but we also need to handle the case where we have
    #  multiple directories. Perhaps if we need to rearrange directory of
    #  directories so that it's above this main. and then we add a param to
    #  supply the path to the sample file, which can be inside the directory
    #  or external. We could also have a switch to say, file is in the
    #  directory.
    logging.debug('Sequence Processing Pipeline Main() Started')

    pipeline = KnightLabPipeline(input_directory,
                                 output_directory,
                                 younger_than=younger_than,
                                 older_than=older_than,
                                 should_filter=should_filter,
                                 nprocs=nprocs)
    pipeline.process_one()

    logging.debug('Sequence Processing Pipeline Main() Completed')

if __name__ == '__main__':
    if True:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.DEBUG,
                            filename='sequence_processing_pipeline.log')

    if sys.argv[3]:
        should_filter = False

    main(sys.argv[1], sys.argv[2], 48, 24, 16, should_filter)
