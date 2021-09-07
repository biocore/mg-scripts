import logging
import shutil
import os
from sequence_processing_pipeline.util import system_call
from time import sleep
from datetime import datetime
from sequence_processing_pipeline.Job import Job


class BCLConvertJob(Job):
    def __init__(self, root_dir, should_filter):
        logging.debug("BCLConvertJob Constructor called")
        self.root_dir = root_dir
        self.should_filter = should_filter
        super().__init__()

    def _get_date(self):
        '''
        Simulates the output of the Linux 'date' command.
        :return:
        '''
        # note the format of timestamp will be similar to:
        # 2021-09-01 21:17:19.005944
        # this is different than the format for 'date':
        # Wed Sep  1 21:11:17 PDT 2021
        # if this is an issue for consumers of this data, we can change it.
        timestamp = str(datetime.now())
        return timestamp

    def _was_job_successful(self, job_num, stdout):
        # process the stdout of sacct into something easy to work with
        results = self._process_sacct_output(job_num, stdout)
        # if one or more component jobs exited with a return code other than
        # 0, return False. Iterate through all jobs rather than return early
        # to write all jobs to logging.
        fail_flag = False
        for job in results:
            # the final element in each job list is the return code
            if job[-1] == '0:0':
                logging.debug("job %s exited with %s" % (str(job), job[-1]))
            else:
                logging.error("job %s exited with %s" % (str(job), job[-1]))
                fail_flag = True

        return fail_flag

    def _process_sacct_output(self, job_num, stdout):
        # treat stdout as a multiline string, and separate into a list of
        # strings.
        l = stdout.split('\n')
        # remove any leading or trailing whitespace in each line
        l = [x.strip() for x in l]
        # remove any line that doesn't begin with the job number.
        # there may be more than one row in the output:
        # 67767+0...
        # 67767+1...
        l = [x for x in l if x.startswith(job_num)]
        # convert each line into a list, so that we have a list of lists.
        # splitting on ' ' will create multiple empty strings due to multiple
        # spaces in between any two columns. Hence, use filter() to filter out
        # those empty strings using 'None' as a parameter.
        l = [list(filter(None, x.split(' '))) for x in l]

        return l

    def run(self, seq_dir, csvfile, base_mask, exp_name, bclconvert_template, polling_wait=300):
        # assume this directory has already been created. we know it has been by this point.

        # I believe root_sequence will always be seq_dir post porting.
        root_sequence = seq_dir

        fastq_output = os.path.join(self.root_dir, 'Data', 'Fastq')
        if self.should_filter == False:
            cmd = ['sbatch',
                   '--parsable',
                   # presumably seq_proc is a defined quality of service level in the system
                   '--qos=seq_proc',
                   '--export=seqdir="%s",outputdir="%s",csvfile="%s",base_mask="%s"' % (
                       seq_dir, fastq_output, csvfile, base_mask),
                   '--job-name=%s' % exp_name,
                   '--partition=long %s' % bclconvert_template
                   ]

            stdout, stderr, return_code = system_call(cmd)
            # there may need some massaging to the output to turn it into a proper job_num
            job_num = stdout  # job_num should remain a string
            job_num_file = os.path.join(fastq_output, job_num) + '.txt'
            with open(job_num_file, 'w') as f:
                f.write('start date == %s\n' % self._get_date())
                f.write('experiment name == %s\n' % exp_name)
                f.write('job_id == %s\n' % job_num)
                f.write(
                    'submit args == sbatch --parsable --qos=seq_proc --export=seqdir="%s",outputdir="%s",csvfile="%s",base_mask="%s" --job-name=%s --partition=long %s' % (
                        seq_dir, fastq_output, csvfile, base_mask, exp_name, bclconvert_template))
            cmd = ['scontrol', 'show', 'job', job_num]
            while True:
                stdout, stderr, return_code = system_call(cmd)
                if return_code == -1:
                    break
                else:
                    logging.debug(
                        "sbatch job %s still in progress. Sleeping %d seconds..." % (job_num, polling_wait))
                    sleep(polling_wait)
            with open(job_num_file, 'a') as f:
                f.write('end date == %s\n' % self._get_date())

            cmd = ['saact', '-X', '-P', '-j', job_num]
            stdout, stderr, return_code = system_call(cmd)
            if self._was_job_successful(job_num, stdout):
                logging.info("Initial processing complete for project %s %s" % (exp_name, root_sequence))
                with open(job_num_file, 'a') as f:
                    # we could also do this for != 0 as well
                    f.write('exit status == 0 for %s\n' % job_num)
                # this may freak out if the file is already there - we need to handle that.
                # we should wrap this for known exceptions.
                shutil.copy(csvfile, fastq_output)
            else:
                logging.info("Initial processing issues %s %s %s" % (exp_name, root_sequence, csvfile))

    # TODO: Not sure yet if we should process processed and alockfile files here or
    #  at the end of process().
    #  touch ${seqdir}/processed && sleep 2
    #  rm ${seqdir}/alockfile
