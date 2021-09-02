from sequence_processing_pipeline.Pipeline import Pipeline


class CMILabPipeline(Pipeline):
    def __init__(self, scan_dir, lab_output_dir, final_output_dir, younger_than=48, older_than=24, should_filter=False, nprocs=16):
        super().__init__(scan_dir, lab_output_dir, final_output_dir, younger_than, older_than, should_filter=should_filter, nprocs=nprocs)
