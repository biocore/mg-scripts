import logging
import sys
from sequence_processing_pipeline.KnightLabPipeline import KnightLabPipeline


def main(seqpath, labname, should_filter, target, file_extension, nprocs, scan_dir, output_dir, younger_than, older_than):
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
                                 younger_than=younger_than,
                                 older_than=older_than,
                                 should_filter=should_filter)
    pipeline.process()

    logging.debug('Sequence Processing Pipeline Main() Completed')

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG) #, filename='sequence_processing_pipeline.log')
    # if should_filter isn't specified as True on the command-line,
    # should_filter should default to false.
    should_filter = False
    main('a', 'b', should_filter, 'd', 'e', 'f', sys.argv[1], 'h')



