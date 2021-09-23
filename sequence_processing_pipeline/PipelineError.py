class PipelineError(Exception):
    def __init__(self, message=None, mailing_list=None):
        self.message = message
        self.mailing_list = mailing_list
        super().__init__(self.message)
