from sequence_processing_pipeline.Pipeline import Pipeline


class KnightLabPipeline(Pipeline):
    def __init__(self, scan_dir, threshold_in_hours=24):
        super().__init__(scan_dir, threshold_in_hours)
