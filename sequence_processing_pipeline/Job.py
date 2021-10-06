from os import makedirs, walk
from os.path import exists, split, join
from sequence_processing_pipeline.PipelineError import PipelineError
from subprocess import Popen, PIPE
from time import sleep
import xml.dom.minidom
import logging


class Job:
    def __init__(self, run_dir, job_name, executable_paths,
                 modules_to_load=None):
        '''
        Base-class to implement Jobs from.
        :param run_dir: The path to a lab run_dir.
        :param job_name: A name for the job. Used to create log files.
        :param executable_paths: A list of executables to validate.
        :param modules_to_load: A list of modules to load before validation.
        '''
        self.run_dir = run_dir
        self.job_name = job_name
        self.stdout_log_path = 'localhost:' + join(self.run_dir,
                                                   f'{self.job_name}.out.log')
        self.stderr_log_path = 'localhost:' + join(self.run_dir,
                                                   f'{self.job_name}.err.log')

        self.script_count = 0

        # For each executable in the list, get its filename and use _which()
        # to see if it can be found. Directly pass an optional list of modules
        # to load before-hand, so that the binary can be found.
        # If the executable can't be found or doesn't have the same path as
        # the version given, raise a PipelineError.
        for executable_path in executable_paths:
            file_path, file_name = split(executable_path)
            # No need to test results. _which() will raise a PipelineError if
            # file_name is a path and the path found does not match. It will
            # also raise a PipelineError if the file could not be found.
            self._which(file_name, modules_to_load=modules_to_load)

    def generate_job_script_path(self):
        '''
        Generate unique paths and filenames to create a Torque job-script with.
        :return: A triplet of strings for job_script_path, out, and err logs.
        '''
        # Note: Not all Job sub-classes will submit jobs to Torque, and not
        # all sub-classes will submit just one job to Torque. This method can
        # be used to generate as many unique paths as is needed for your Job,
        # following a standard convention.
        self.script_count += 1

        job_script_path = join(self.run_dir,
                               f'{self.job_name}_{self.script_count}.sh')

        output_log_path = join('localhost:',
                               self.run_dir,
                               f'{self.job_name}_{self.script_count}.out.log')

        error_log_path = join('localhost:',
                              self.run_dir,
                              f'{self.job_name}_{self.script_count}.err.log')

        return (job_script_path, output_log_path, error_log_path)

    def run(self):
        '''
        Since a Job object can encapsulate one or more qsub() or system()
        calls, the base run() method remains unimplemented. It is the job of
        each sub-class to define what needs to be run in order to generate the
        expected output.
        '''
        raise PipelineError("Base class run() method not implemented.")

    def _which(self, file_path, modules_to_load=None):
        '''
        Returns file_path if file_path exists and file_path is a full path.
        Otherwise returns the path to a file named 'file_path' found in PATH.
        :param file_path: The path of the executable to find.
        :param modules_to_load: A list of Linux module names to load.
        :return: A path to 'file_name'.
        '''
        tmp = split(file_path)
        # remove any elements that are empty string.
        tmp = [x for x in tmp if x]

        isPath = True if len(tmp) > 1 else False

        cmd = 'which ' + file_path

        if modules_to_load:
            cmd = 'module load ' + ' '.join(modules_to_load) + ';' + cmd

        results = self._system_call(cmd)
        result = results['stdout'].strip()

        if not result:
            raise PipelineError("File '%s' does not exist." % file_path)

        if isPath is True and result != file_path:
            raise PipelineError(f"Found path '{result} does not match "
                                f"{file_path}")

        return result

    def _file_check(self, file_path):
        if exists(file_path):
            logging.debug("file '%s' exists." % file_path)
            return True
        else:
            raise PipelineError("file '%s' does not exist." % file_path)

    def _find_files(self, search_path):
        lst = []
        for root, dirs, files in walk(search_path):
            lst += [join(root, x) for x in files]
        return lst

    def _directory_check(self, directory_path, create=False):
        if exists(directory_path):
            logging.debug("directory '%s' exists." % directory_path)
        else:
            if create:
                try:
                    makedirs(directory_path, exist_ok=True)
                except OSError as e:
                    # this is a known potential error. Re-raise it as a
                    # PipelineError, so it gets handled in the same location
                    # as the others.
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
        '''
        Submit a Torque job script and optionally wait for it to finish.
        :param script_path: The path to a Torque job (bash) script.
        :param qsub_parameters: Optional parameters for qsub.
        :param script_parameters: Optional parameters for your job script.
        :param wait: Set to False to submit job and not wait.
        :return: Dictionary containing the job's id, name, status, and
        elapsed time. Raises PipelineError if job could not be submitted or
        if job was unsuccessful.
        '''
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
        results = self._system_call(cmd)
        stdout = results['stdout']

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

            results = self._system_call("qstat -x %s" % job_id,
                                        allow_return_codes=[153])
            stdout = results['stdout']
            stderr = results['stderr']

            logging.debug("QSTAT STDOUT: [%s]\n" % stdout)
            logging.debug("QSTAT STDERR: [%s]\n" % stderr)

            if stderr.startswith("qstat: Unknown Job Id Error"):
                break

            # update job_info
            # even though job_id doesn't change, use this update as evidence
            # job_id was once in the queue.
            doc = xml.dom.minidom.parseString(stdout)
            jobs = doc.getElementsByTagName("Job")
            for some_job in jobs:
                tmp = some_job.getElementsByTagName("Job_Id")[0]
                some_job_id = tmp.firstChild.nodeValue
                if some_job_id == job_id:
                    job_info['job_id'] = some_job_id
                    tmp = some_job.getElementsByTagName("Job_Name")[0]
                    job_info['job_name'] = tmp.firstChild.nodeValue
                    tmp = some_job.getElementsByTagName("job_state")[0]
                    job_info['job_state'] = tmp.firstChild.nodeValue
                    # we cannot use exit_status, as it may appear for regular
                    # jobs but apparently it doesn't appear for job arrays.
                    # Use status = 'C' instead.
                    tmp = some_job.getElementsByTagName("etime")[0]
                    job_info['elapsed_time'] = tmp.firstChild.nodeValue
                    break

            logging.debug("Job info: %s" % job_info)

            # if job is completed after having run or exited after having
            # run, then stop waiting.
            if job_info['job_state'] in ['C', 'E']:
                break

            # check once every minute - job info records should stay in the
            # queue that long after finishing.
            sleep(30)

        if job_info['job_id']:
            # job was once in the queue
            if job_info['job_state'] == 'C':
                # job completed successfully
                return job_info
            else:
                # job exited unsuccessfully
                status = job_info['job_state']
                raise PipelineError("job %s exited with status %s" % (job_id,
                                                                      status))
        else:
            # job was never in the queue - return an error.
            raise PipelineError("job %s never appeared in the queue." % job_id)
