from sequence_processing_pipeline.Pipeline import Pipeline


class KnightLabPipeline(Pipeline):
    def __init__(self, scan_dir, younger_than=48, older_than=24, should_filter=False):
        super().__init__(scan_dir, younger_than, older_than, should_filter=should_filter)
