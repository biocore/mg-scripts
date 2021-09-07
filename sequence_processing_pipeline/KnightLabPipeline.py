from sequence_processing_pipeline.Pipeline import Pipeline


class KnightLabPipeline(Pipeline):
    def __init__(self, scan_dir, lab_output_dir, younger_than=48,
                 older_than=24, should_filter=False, nprocs=16):
        # for labs that don't define a final_output_dir, set it equal to
        # lab_output_dir, for now.
        super().__init__(scan_dir, lab_output_dir, lab_output_dir,
                         younger_than, older_than, should_filter=should_filter,
                         nprocs=nprocs)
