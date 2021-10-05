import logging
from sequence_processing_pipeline.Pipeline import Pipeline


def main(configuration_file_path, sample_sheet_path, input_directory,
         output_directory, run_id):
    '''
    Provide an easy entry point for users who would import the package and use
    it. It defines the parameters needed to process a directory of data.
    :param sample_sheet_path: A path to a sample-sheet for the BCL files.
    :param input_directory: A directory containing BCL files.
    :param output_directory: A directory to store the Pipeline's products.
    '''
    logging.debug('Sequence Processing Pipeline main() Started')

    pipeline = Pipeline(configuration_file_path, input_directory,
                        output_directory, run_id)
    pipeline.process(sample_sheet_path)

    logging.debug('Sequence Processing Pipeline main() Completed')
