import logging
from sequence_processing_pipeline.Job import Job
from time import sleep
from sequence_processing_pipeline.PipelineError import PipelineError


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

    def qsub(self, script_path, qsub_parameters='', script_parameters=''):
        # -w e: verify options and abort if there's an error.
        # we now define queue (-q long) in the job script.
        cmd = 'qsub -w e %s %s %s' % (qsub_parameters,
                                      script_path,
                                      script_parameters)

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
            stdout, stderr = self._system_call("qstat %s" % job_id,
                                               allow_return_codes=[153])

            # assume stdout appears like this:
            # Job ID      Username  Queue         Jobname     Limit  State  Elapsed
            # ------      --------  -----         -------     -----  -----  -------
            # 385         jqpublic  sun-medium    STDIN       06:00  Run    00:03
            # or like this:
            # qstat: Unknown Job Id Error 992223.barnacle.ucsd.edu

            if stdout.startswith("qstat: Unknown Job Id Error"):
                break

            # break stdout into multiple lines, and extract the line that begins
            # with the job_id.
            lines = stdout.strip().split('\n')
            # remove any beginning or trailing whitespace from each line.
            lines = [x.strip() for x in lines]
            # if the line begins with job_id, then that is the one and only one
            # line we want. To find it, we'll use a list comprehension instead
            # of breaking a for loop for speed and convenient syntax.
            # Once we find the line containing job_id, we'll split the line on
            # whitespace. Since this is a list containing exactly one list, we'll
            # remove the outer list by taking the element 0.
            lines = [x.split(' ') for x in lines if x.startswith(job_id)][0]
            # remove any element from the list that's ''.
            lines = [x for x in lines if x]

            # update job_info
            # even though job_id doesn't change, use this update as evidence
            # job_id was once in the queue.
            job_info['job_id'] = job_id
            job_info['job_name'] = lines[3]
            job_info['status'] = lines[5]
            job_info['elapsed_time'] = lines[6]

            logging.debug("Job status for %s" % (job_info))

            if job_info['status'] in ['C', 'E']:
                # if job is in completed or exiting state, then exit.
                # this exits earlier than waiting for the entry to
                # fall out of qstat, and with a known status.
                break

            # check once every minute - job info records should stay in the
            # queue that long after finishing.
            sleep(60)

        if job_info['job_id']:
            # job was once in the queue
            if job_info['status'] == 'C':
                # job completed successfully
                return job_info
            else:
                # job exited unsuccessfully
                raise PipelineError("job %s exited with status %s" % (job_id, job_info['status']))
        else:
            # job was never in the queue - return an error.
            raise PipelineError("job %s never appeared in the queue." % job_id)
