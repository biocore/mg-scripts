class PipelineError(Exception):
    def __init__(self, message=None, mailing_list=None):
        self.message = message
        self.mailing_list = mailing_list
        super().__init__(self.message)


class JobFailedError(PipelineError):
    # Occurs when a successfully-submitted job as failed.
    def __init__(self, message=None, mailing_list=None):
        self.message = message
        self.mailing_list = mailing_list
        super().__init__(self.message)


class ExecFailedError(PipelineError):
    # Occurs when an executed command returns w/an error, which is defined as
    # the command returning a value not zero and not an acceptable non-zero
    # value.
    def __init__(self, message=None, mailing_list=None):
        self.message = message
        self.mailing_list = mailing_list
        super().__init__(self.message)


class LogParsingError(PipelineError):
    # Occurs when an executed command returns w/an error, which is defined as
    # the command returning a value not zero and not an acceptable non-zero
    # value. May or may not be useful.
    def __init__(self, message=None, mailing_list=None):
        self.message = message
        self.mailing_list = mailing_list
        super().__init__(self.message)
