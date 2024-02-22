from itertools import zip_longest
from os import makedirs, walk
from os.path import basename, exists, split, join
from sequence_processing_pipeline.PipelineError import (PipelineError,
                                                        JobFailedError,
                                                        ExecFailedError)
from subprocess import Popen, PIPE
from time import sleep
import logging
from inspect import stack


class Job:
    def __init__(self, root_dir, output_path, job_name, executable_paths,
                 max_array_length, modules_to_load=None):
        """
        Base-class to implement Jobs from.
        :param job_name: A name for the job. Used to create log files.
        :param root_dir: The path to a Job's root input directory.
        :param output_path: The root path to store all job products.
        :param executable_paths: A list of executables to validate.
        :param modules_to_load: A list of modules to load before validation.
        """
        self.job_name = job_name
        self.root_dir = root_dir
        self._directory_check(self.root_dir, create=False)
        self.force_job_fail = False

        self.output_path = join(output_path, self.job_name)
        self._directory_check(self.output_path, create=True)

        self.log_path = join(self.output_path, 'logs')
        self._directory_check(self.log_path, create=True)

        self.modules_to_load = modules_to_load
        self.max_array_length = max_array_length

        self.script_count = 0

        self.bypass_exec_check = ['bcl-convert']

        self.suffix = None

        # checking if this is running as part of the unittest
        # https://stackoverflow.com/a/25025987
        self.is_test = True if [
            x for x in stack() if 'unittest' in x.filename] else False

        # For each executable in the list, get its filename and use _which()
        # to see if it can be found. Directly pass an optional list of modules
        # to load before-hand, so that the binary can be found.
        # If the executable can't be found or doesn't have the same path as
        # the version given, raise a PipelineError.
        for executable_path in executable_paths:
            file_path, file_name = split(executable_path)

            # bcl-convert module is not installed on the node this test gets
            # run on, so forego it entirely.

            # Some modules such as bcl-convert are not available to the
            # servers running this check. However they are still available
            # on Barnacle nodes and this is what's important. For now simply
            # bypass the check for known situations.
            for name in self.bypass_exec_check:
                if name in file_name:
                    continue

            # No need to test results. _which() will raise a PipelineError if
            # file_name is a path and the path found does not match. It will
            # also raise a PipelineError if the file could not be found.
            if not self.is_test:
                self._which(file_name, modules_to_load=self.modules_to_load)

    def run(self):
        """
        Since a Job object can encapsulate one or more submit_job() or system()
        calls, the base run() method remains unimplemented. It is the job of
        each sub-class to define what needs to be run in order to generate the
        expected output.
        """
        raise PipelineError("Base class run() method not implemented.")

    def parse_logs(self):
        raise PipelineError("Base class parse_logs() method not implemented.")

    def _which(self, file_path, modules_to_load=None):
        """
        Returns file_path if file_path exists and file_path is a full path.
        Otherwise returns the path to a file named 'file_path' found in PATH.
        :param file_path: The path of the executable to find.
        :param modules_to_load: A list of Linux module names to load.
        :return: A path to 'file_name'.
        """
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

    def _system_call(self, cmd, allow_return_codes=[], callback=None):
        """
        Call command and return (stdout, stderr, return_value)
        :param cmd: The string containing the command to be run, or a sequence
                    of strings that are the tokens of the command.
        :param allow_return_codes: optional user-defined list of successful
                    return codes in addition to 0.
        :param callback: optional function taking two parameters (id, status)
                         that is called when a running process's status is
                         changed.
        :return: a dictionary containing stdout, stderr, and return_code as
                 key/value pairs.
        """
        proc = Popen(cmd, universal_newlines=True, shell=True,
                     stdout=PIPE, stderr=PIPE)

        if callback is not None:
            callback(jid=proc.pid, status='RUNNING')

        # Communicate pulls all stdout/stderr from the PIPEs
        # This call blocks until the command is done
        stdout, stderr = proc.communicate()
        return_code = proc.returncode

        logging.debug("stdout: %s" % stdout)
        logging.debug("stderr: %s" % stderr)
        logging.debug("return code: %s" % return_code)

        acceptable_return_codes = [0] + allow_return_codes

        if return_code not in acceptable_return_codes:
            if callback is not None:
                callback(jid=proc.pid, status='ERROR')
            msg = (
                'Execute command-line statement failure:\n'
                f'Command: {cmd}\n'
                f'return code: {return_code}\n'
                f'stdout: {stdout}\n'
                f'stderr: {stderr}\n')
            logging.error(msg)
            raise ExecFailedError(message=msg)

        if callback is not None:
            callback(jid=proc.pid, status='COMPLETED')

        return {'stdout': stdout, 'stderr': stderr, 'return_code': return_code}

    def submit_job(self, script_path, job_parameters=None,
                   script_parameters=None, wait=True,
                   exec_from=None, callback=None):
        """
        Submit a Torque job script and optionally wait for it to finish.
        :param script_path: The path to a Torque job (bash) script.
        :param job_parameters: Optional parameters for scheduler submission.
        :param script_parameters: Optional parameters for your job script.
        :param wait: Set to False to submit job and not wait.
        :param exec_from: Set working directory to execute command from.
        :param callback: Set callback function that receives status updates.
        :return: Dictionary containing the job's id, name, status, and
        elapsed time. Raises PipelineError if job could not be submitted or
        if job was unsuccessful.
        """
        if job_parameters:
            cmd = 'sbatch %s %s' % (job_parameters, script_path)
        else:
            cmd = 'sbatch %s' % (script_path)

        if script_parameters:
            cmd += ' %s' % script_parameters

        if exec_from:
            cmd = f'cd {exec_from};' + cmd

        logging.debug("job scheduler call: %s" % cmd)

        if self.force_job_fail:
            raise JobFailedError("This job died.")

        # if system_call does not raise a PipelineError(), then the scheduler
        # successfully submitted the job. In this case, it should return
        # the id of the job in stdout.
        results = self._system_call(cmd)
        stdout = results['stdout']

        job_id = stdout.strip().split()[-1]

        job_info = {'job_id': None, 'job_name': None, 'job_state': None,
                    'elapsed_time': None}
        # Just to give sometime for everything to be set up properly
        sleep(5)
        while wait:
            result = self._system_call(f"sacct -P -n --job {job_id} --format "
                                       "JobID,JobName,State,Elapsed,ExitCode")

            if result['return_code'] != 0:
                # sacct did not successfully submit the job.
                raise ExecFailedError(result['stderr'])

            # [-1] remove the extra \n
            jobs_data = result['stdout'].split('\n')[:-1]
            states = dict()
            estatuses = dict()
            for i, jd in enumerate(jobs_data):
                jid, jname, jstate, etime, estatus = jd.split('|')
                if jid.endswith('.extern') or jid.endswith('.batch'):
                    continue

                if i == 0:
                    job_info['job_id'] = jid
                    job_info['job_name'] = jname
                    job_info['elapsed_time'] = etime
                    job_info['exit_status'] = estatus

                if jstate not in states:
                    states[jstate] = 0
                states[jstate] += 1

                if estatus not in estatuses:
                    estatuses[estatus] = 0
                estatuses[estatus] += 1

            job_info['job_state'] = f'{states}'
            job_info['exit_status'] = f'{estatuses}'

            if callback is not None:
                callback(jid=job_id, status=f'{states}')

            logging.debug("Job info: %s" % job_info)

            # if job is completed after having run or exited after having
            # run, then stop waiting.
            if not set(states) - {'COMPLETED', 'FAILED', 'CANCELLED'}:
                break

            sleep(5)

        if job_info['job_id'] is not None:
            # job was once in the queue
            if callback is not None:
                callback(jid=job_id, status=job_info['job_state'])

            if set(states) == {'COMPLETED'}:
                if 'exit_status' in job_info:
                    if set(estatuses) == {'0:0'}:
                        # job completed successfully
                        return job_info
                    else:
                        exit_status = job_info['exit_status']
                        raise JobFailedError(f"job {job_id} exited with exit_"
                                             f"status {exit_status}")
                else:
                    # with no other info, assume job completed successfully
                    return job_info
            else:
                # job exited unsuccessfully
                raise JobFailedError(f"job {job_id} exited with status "
                                     f"{job_info['job_state']}")
        else:
            # job was never in the queue - return an error.
            if callback is not None:
                callback(jid=job_id, status='ERROR')

            raise JobFailedError(f"job {job_id} never appeared in the "
                                 "queue.")

    def _group_commands(self, cmds):
        # break list of commands into chunks of max_array_length (Typically
        # 1000 for Torque job arrays). To ensure job arrays are never more
        # than 1000 jobs long, we'll chain additional commands together, and
        # evenly distribute them amongst the first 1000.
        cmds.sort()
        chunks = [cmds[i:i + self.max_array_length] for i in
                  range(0, len(cmds), self.max_array_length)]

        results = []

        # create a chained command by taking one command from each list.
        # zip_longest() allows us to handle lists of different lengths, as the
        # last chunk will always be of different length than 1000.
        for tuple in zip_longest(*chunks):
            # zip_longest() pads shorter lists with None. In our case, we
            # don't want an additional command named 'None'.
            chained_cmd = [x for x in list(tuple) if x is not None]
            chained_cmd = ';'.join(chained_cmd)
            results.append(chained_cmd)

        return results

    # assume for now that a corresponding zip file exists for each html
    # file found. Assume for now that all html files will be found in a
    # 'filtered_sequences' or 'trimmed_sequences' subdirectory.
    #
    # verify that the entire list of sample-ids found match what's
    # expected. Since the list of expected ids is very small, we'll
    # perform an exact comparison.

    def audit(self, sample_ids):
        """
        Audit the results of a run.
        :param sample_ids: A list of sample-ids that require results.
        :return: A list of sample-ids that were not found.
        """
        files_found = []

        if self.suffix is None:
            raise PipelineError("Audit() method called on base Job object.")

        for root, dirs, files in walk(self.output_path):
            if 'zero_files' in root:
                continue
            files_found += [join(root, x) for x in files if
                            x.endswith(self.suffix)]

        found = []
        for sample_id in sample_ids:
            for found_file in files_found:
                # the trailing underscore is important as it can be assumed
                # that all fastq.gz files will begin with sample_id followed
                # by an '_', and then one or more additional parameters
                # separated by underscores. This substring is unlikely to be
                if basename(found_file).startswith('%s_' % sample_id):
                    found.append(sample_id)
                    break

        return sorted(list(set(found) ^ set(sample_ids)))

    def _toggle_force_job_fail(self):
        if self.force_job_fail is True:
            self.force_job_fail = False
        else:
            self.force_job_fail = True

        return self.force_job_fail
