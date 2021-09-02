import logging
import sys
from sequence_processing_pipeline.KnightLabPipeline import KnightLabPipeline


def main(should_filter, nprocs, scan_dir, lab_output_dir, younger_than, older_than):
    # the purpose of main() is to provide an easy entry point for users who
    # would import the package and use it. It defines the parameters needed
    # to process a directory of data.
    #
    # The real work is done in the Pipeline object. We need to clean this
    # up so fatal errors and what-not are handled in the right places.
    #
    # labname can be fed as a parameter to Pipeline, or it can be used to
    # select a subclass of Pipeline tailored to each lab.
    logging.debug('Sequence Processing Pipeline Main() Started')

    pipeline = KnightLabPipeline(scan_dir,
                                 lab_output_dir,
                                 younger_than=younger_than,
                                 older_than=older_than,
                                 should_filter=should_filter,
                                 nprocs=nprocs)
    pipeline.process()

    logging.debug('Sequence Processing Pipeline Main() Completed')

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG) #, filename='sequence_processing_pipeline.log')
    # if should_filter isn't specified as True on the command-line,
    # should_filter should default to false.
    should_filter = False
    main(should_filter, 16, sys.argv[1], sys.argv[2], 48, 24)




