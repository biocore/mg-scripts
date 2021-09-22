import logging
from subprocess import Popen, PIPE
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import exists
from os import makedirs


class Job:
    def __init__(self):
        logging.debug("Job Constructor called.")

    def _directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PinelineError, so it gets handled in the same location as the
                    # others.
                    raise PipelineError(str(e))
            else:
                raise PipelineError("directory_path '%s' does not exist." % directory_path)

    def _system_call(self, cmd, allow_return_codes=[]):
        """Call command and return (stdout, stderr, return_value)

        Parameters
        ----------
        cmd : str or iterator of str
            The string containing the command to be run, or a sequence of strings
            that are the tokens of the command.

        Returns
        -------
        str, str, int
            - The standard output of the command
            - The standard error of the command
            - The exit status of the command

        Notes
        -----
        This function is ported from QIIME (http://www.qiime.org), previously named
        qiime_system_call. QIIME is a GPL project, but we obtained permission from
        the authors of this function to port it to Qiita and keep it under BSD
        license.
        """
        logging.debug("Job _system_call() method called.")

        proc = Popen(cmd, universal_newlines=True, shell=True, stdout=PIPE, stderr=PIPE)
        # Communicate pulls all stdout/stderr from the PIPEs
        # This call blocks until the command is done
        stdout, stderr = proc.communicate()
        return_code = proc.returncode

        logging.debug("stdout: %s" % stdout)
        logging.debug("stderr: %s" % stderr)
        logging.debug("return code: %s" % return_code)

        acceptable_return_codes = [0] + allow_return_codes

        if return_code not in acceptable_return_codes:
            s = "Execute command-line statement failure:\n"
            s += "return code: %s\n" % return_code
            s += "stdout: %s\n" % stdout
            s += "stderr: %s\n" % stderr
            logging.error(s)
            raise PipelineError(message=s)

        return {'stdout': stdout, 'stderr': stderr, 'return_code': return_code}


class lsJob(Job):
    def __init__(self, path):
        super().__init__()
        logging.debug("lsJob Constructor called.")
        self.path = path

    def run(self):
        logging.debug("lsJob run() method called.")

        cmd = 'ls %s' % self.path
        logging.debug(cmd)
        return self._system_call(cmd)
