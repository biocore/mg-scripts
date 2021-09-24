from os import makedirs
from os.path import exists
from re import search
from sequence_processing_pipeline.PipelineError import PipelineError
from subprocess import Popen, PIPE
from time import sleep
import logging


class Job:
    def __init__(self):
        pass

    def run(self):
        raise PipelineError("Base class run() method not implemented.")

    def _directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PipelineError, so it gets handled in the same location as
                    # the others.
                    raise PipelineError(str(e))
            else:
                raise PipelineError(
                    "directory_path '%s' does not exist." % directory_path)

    def _system_call(self, cmd, allow_return_codes=[]):
        """Call command and return (stdout, stderr, return_value)

        Parameters
        ----------
        cmd : str or iterator of str
            The string containing the command to be run, or a sequence of
            strings that are the tokens of the command.

        Returns
        -------
        str, str, int
            - The standard output of the command
            - The standard error of the command
            - The exit status of the command

        Notes
        -----
        This function is ported from QIIME (http://www.qiime.org), previously
        named qiime_system_call. QIIME is a GPL project, but we obtained
        permission from the authors of this function to port it to Qiita and
        keep it under BSD license.
        """
        logging.debug("Job _system_call() method called.")

        proc = Popen(cmd, universal_newlines=True, shell=True,
                     stdout=PIPE, stderr=PIPE)
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

    def qsub(self, script_path, qsub_parameters=None,
             script_parameters=None, wait=True):
        # -w e: verify options and abort if there's an error.
        # we now define queue (-q long) in the job script.
        # cmd = 'qsub -w e %s %s %s' % (qsub_parameters,
        if qsub_parameters:
            cmd = 'qsub %s %s' % (qsub_parameters, script_path)
        else:
            cmd = 'qsub %s' % (script_path)

        if script_parameters:
            cmd += ' %s' % script_parameters

        logging.debug("QSUB call: %s" % cmd)
        # if system_call does not raise a PipelineError(), then the qsub
        # successfully submitted the job. In this case, qsub should return
        # the id of the job in stdout.
        stdout, stderr = self._system_call(cmd)

        # extract the first line only. This will be the fully-qualified Torque
        # job id.
        job_id = stdout.split('\n')[0]

        # job_info acts as a sentinel value, as well as a return value.
        # if job_id is not None, then that means the job has appeared in qstat
        # at least once. This helps us distinguish between two states -
        # one where the loop is starting and the job has never appeared in
        # qstat yet, and the second where we slept too long in between
        # checking and the job has disappeared from qstat().
        job_info = {'job_id': None, 'job_name': None, 'status': None,
                    'elapsed_time': None}

        while wait:
            # wait for the job to complete. Periodically check on the status
            # of the job using the 'qstat' command before sleep()ing again.
            # Recently completed and exited jobs appear in qstat for only a
            # short period of time. If job_id is not found, treat this as a\
            # job finishing, rather than an error. qstat returns code 153 when
            # this occurs, so we'll allow it in this instance.

            stdout, stderr = self._system_call("qstat -x %s" % job_id,
                                               allow_return_codes=[153])

            logging.debug("QSTAT STDOUT: %s\n" % stdout)
            logging.debug("QSTAT STDERR: %s\n" % stderr)

            if stdout.startswith("qstat: Unknown Job Id Error"):
                break

            # update job_info
            # even though job_id doesn't change, use this update as evidence
            # job_id was once in the queue.
            job_info['job_id'] = search(
                r"<Job_Id>(.*?)</Job_Id>", stdout).group(1)
            job_info['job_name'] = search(
                r"<Job_Name>(.*?)</Job_Name>", stdout).group(1)
            job_info['status'] = search(
                r"<job_state>(.*?)</job_state>", stdout).group(1)
            job_info['elapsed_time'] = search(
                r"<etime>(.*?)</etime>", stdout).group(1)
            exit_status = search(
                r"<exit_status>(.*?)</exit_status>", stdout)
            if exit_status:
                job_info['exit_status'] = exit_status.group(1)

            logging.debug("Job info: %s" % job_info)

            if 'exit_status' in job_info:
                break

            # check once every minute - job info records should stay in the
            # queue that long after finishing.
            sleep(30)

        if job_info['job_id']:
            # job was once in the queue
            if job_info['exit_status'] == '0':
                # job completed successfully
                return job_info
            else:
                # job exited unsuccessfully
                status = job_info['exit_status']
                raise PipelineError("job %s exited with status %s" % (job_id,
                                                                      status))
        else:
            # job was never in the queue - return an error.
            raise PipelineError("job %s never appeared in the queue." % job_id)
