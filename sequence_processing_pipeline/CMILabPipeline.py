from sequence_processing_pipeline.Pipeline import Pipeline


class CMILabPipeline(Pipeline):
    def __init__(self, scan_dir, younger_than=24):
        super().__init__(scan_dir, younger_than)
