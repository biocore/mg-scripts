from sequence_processing_pipeline.Pipeline import Pipeline


class CMILabPipeline(Pipeline):
    def __init__(self, input_directory, output_directory,
                 final_output_directory, younger_than=48,
                 older_than=24, nprocs=16):
        super().__init__(input_directory, output_directory, final_output_directory,
                         younger_than, older_than, nprocs=nprocs)
