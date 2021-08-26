class PipelineError(Exception):
    def __init__(self, mailing_list, message=None):
        self.mailing_list = mailing_list
        self.message = message
        super().__init__(self.message)

