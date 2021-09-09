import logging
from sequence_processing_pipeline.util import system_call
from time import sleep


class Job:
    def __init__(self, polling_wait=300):
        logging.debug("Base Job Constructor called")
        self.polling_wait = polling_wait

    def execute_sbatch_job_and_wait(self, cmd):
        stdout, stderr, return_code = system_call(cmd)
        job_id = stdout
        cmd = ['scontrol', 'show', 'job', job_id]
        while True:
            stdout, stderr, return_code = system_call(cmd)
            if return_code == -1:
                break
            else:
                logging.debug("sbatch job %s still in progress. Sleeping %d seconds..." % (job_id, self.polling_wait))
                sleep(self.polling_wait)

        return job_id
