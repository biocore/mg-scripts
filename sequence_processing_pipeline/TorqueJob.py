import logging
from sequence_processing_pipeline.Job import Job
from time import sleep
from sequence_processing_pipeline.PipelineError import PipelineError
import re


class TorqueJob(Job):
    def __init__(self):
        super().__init__()
        logging.debug("TorqueJob Constructor called.")
        logging.debug("TorqueJob Constructor exiting.")

    def run(self):
        logging.debug("TorqueJob run() method called.")
        self.qsub('myqueue', 'nooptions', '/mypath', 'noparameters')
        logging.debug("TorqueJob run() method exiting.")

    def qsub2(self):
        options = '-N %s -l h_vmem=4G -pe smp <num_slots> -o outputlogfile -e errorlogfile'

    def qsub(self, script_path, qsub_parameters=None, script_parameters=None):
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

        # extract the first line only, and extract only the first component
        # of the line. This will be the Torque job id.
        job_id = stdout.split('\n')[0].split('.')[0]

        job_info = {'job_id': None, 'job_name': None, 'status': None,
                    'elapsed_time': None}

        while True:
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

            # TODO: Use XML parsing package?
            job_id = re.search(r"<Job_Id>(.*?)</Job_Id>", stdout).group(1)
            job_name = re.search(r"<Job_Name>(.*?)</Job_Name>", stdout).group(1)
            status = re.search(r"<job_state>(.*?)</job_state>", stdout).group(1)
            start_time = re.search(r"<ctime>(.*?)</ctime>", stdout).group(1)
            elapsed_time = re.search(r"<etime>(.*?)</etime>", stdout).group(1)
            exit_status = re.search(r"<exit_status>(.*?)</exit_status>", stdout)

            # update job_info
            # even though job_id doesn't change, use this update as evidence
            # job_id was once in the queue.
            job_info['job_id'] = job_id
            job_info['job_name'] = job_name
            job_info['status'] = status
            job_info['elapsed_time'] = elapsed_time
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
                raise PipelineError("job %s exited with status %s" % (job_id, job_info['exit_status']))
        else:
            # job was never in the queue - return an error.
            raise PipelineError("job %s never appeared in the queue." % job_id)
