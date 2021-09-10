import logging
from subprocess import Popen, PIPE
from sequence_processing_pipeline.PipelineError import PipelineError


class Job:
    def __init__(self):
        logging.debug("Job Constructor called.")
        logging.debug("Job Constructor exiting.")

    def run(self):
        logging.debug("Job run() method called.")
        logging.debug("Job run() method exiting.")

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
        proc = Popen(cmd, universal_newlines=True, shell=True, stdout=PIPE, stderr=PIPE)
        # Communicate pulls all stdout/stderr from the PIPEs
        # This call blocks until the command is done
        stdout, stderr = proc.communicate()
        return_code = proc.returncode

        acceptable_return_codes = [0] + allow_return_codes

        if return_code not in acceptable_return_codes:
            s = "Execute command-line statement failure:\n"
            s += "return code: %s\n" % return_code
            s += "stdout: %s\n" % stdout
            s += "stderr: %s\n" % stderr
            logging.error(s)
            raise PipelineError(message=s)

        return stdout, stderr


class lsJob(Job):
    def __init__(self, path):
        super().__init__()
        logging.debug("lsJob Constructor called.")
        self.path = path
        logging.debug("lsJob Constructor exiting.")

    def run(self):
        logging.debug("lsJob run() method called.")
        out, err, rc = self._system_call(['ls', self.path])
        logging.debug("ls returned: %d" % rc)
        logging.debug("ls stdout: %s" % out)
        logging.debug("ls stderr: %s" % err)
        logging.debug("lsJob run() method exiting.")


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    job = Job()
    job.run()

    logging.debug("Hello World!")

    job = lsJob('/tmp')
    job.run()
