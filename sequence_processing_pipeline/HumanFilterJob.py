import logging
import os
from sequence_processing_pipeline.Job import Job


class HumanFilterJob(Job):
    def __init__(self, SequenceDirectoryObject, nprocs, job_owner_home, email_list, chemistry, output_dir, final_output_dir):
        self.sdo = SequenceDirectoryObject
        self.seq_dir = self.sdo.seq_dir
        self.csv_files = []
        self.nprocs = nprocs
        self.home = job_owner_home
        self.email_list = email_list
        self.chemistry = chemistry
        self.output_dir = output_dir
        self.final_output_dir = final_output_dir
        super().__init__()

    def _get_run_config_info(self, sample_sheet_path):
        with open(sample_sheet_path, 'r') as f:
            # we can (and should) open up this file as a proper INI
            # file. However it may have trailing ',' characters,
            # which may impede this. Hence, open it as a regular
            # file and clean it up before processing.
            lines = f.readlines()
            # first, let's strip out \n and \r\n and any leading or
            # trailing whitespaces
            lines = [x.strip() for x in lines]
            # second, let's strip out trailing ',' characters, if
            # they're present. We'll keep any leading ',' characters
            # or in-line ',' characters.
            lines = [x.rstrip(',') for x in lines]
            # lastly, there have been cases where lines contain only
            # ',' characters. Hence, look for and remove these now
            # empty lines before continuing.
            lines = [x for x in lines if x]
            # since the file is already in memory, we won't use INI
            # library to parse it, (for now). Although it would be
            # cleaner.
            metadata = []
            sentinel = False
            for i in range(0, len(lines)):
                if lines[i] == '[Bioinformatics]':
                    # when Bioinformatics section is found, start
                    # copying lines to the list buffer, but don't
                    # copy the header itself.
                    sentinel = True
                elif lines[i].startswith('['):
                    # when the next header is found, stop copying
                    # lines to the list buffer. Don't include this
                    # header line, either.
                    sentinel = False
                elif sentinel == True:
                    # this should be a line in between
                    # [Bioinformatics] and the next section. Copy it
                    # to the buffer.
                    metadata.append(lines[i])
                # if the end of the file is reached before the next
                # header is found, that means there were no more
                # headers and that's okay too.

            # remove duplicate lines (this appears to be an issue in
            # the original bash scripts.) and sort.
            l = list(set(metadata)).sort()

            # before returning, split each line into a list, so it's
            # easier to iterate through each line of values.
            return [x.split(',') for x in l]

    def _split_list(self, original_list, lines_per_file, file_prefix, file_extension, destination_path):
        def _chunk_list(some_list, n):
            # taken from https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
            for i in range(0, len(some_list), n):
                yield some_list[i:i + n]

        new_files = []
        count = 0
        for chunk in _chunk_list(original_list, lines_per_file):
            # create the filename for the new split file
            s = '%s%d.%s' % (file_prefix, count, file_extension)
            # prepend the path of the original file to the new split files
            s = os.path.join(destination_path, s)
            with open(s, 'w') as f2:
                for line in chunk:
                    f2.write(line)
            new_files.append(s)
            count += 1

        return new_files

    def create_qsub_file(self, seq_dir, project, file_base, home_dir,
                            email_list, trim_file, slurm_array_task_id,
                            chemistry, output_dir, adapter_a, adapter_A,
                            a_trim, h_filter, qiita_proj, final_output_dir):
        lines = []
        lines.append("#!/bin/bash")
        lines.append("#SBATCH --job-name={}_%A_%a".format(project))
        lines.append("#SBATCH --ntasks=%d" % self.nprocs)
        lines.append("#SBATCH --ntasks-per-node=%d" % self.nprocs)
        lines.append("#SBATCH --export=ALL")
        # assuming for now that 72 hours for all jobs is appropriate
        lines.append("#SBATCH --time=72:00:00")
        # assuming for now that 6G is proper for all jobs.
        lines.append("#SBATCH --mem-per-cpu=6G")
        lines.append("#SBATCH --output={}/filter_jobs/%x_%A_%a.out".format(home_dir))
        lines.append("#SBATCH --error={}/filter_jobs/%x_%A_%a.err".format(home_dir))
        lines.append("#SBATCH --mail-type=ALL")
        lines.append("#SBATCH --mail-user=%s" % ','.join(email_list))
        # assume for now this hardcoded path is acceptable
        lines.append("source ~/miniconda3/bin/activate test_env_2")
        # assume for now the backslash below is intentional
        lines.append("file=%s\%s" % (trim_file, slurm_array_task_id))
        lines.append('cmd="sh ~/seq_proc_dev/fpmmp.sh -d %s -D /Data/Fastq -S %s\%s -p %s -C %s -c %d -o %s -O %s -a %s -A %s -g %s -G %s -q %s -f %s"' % (seq_dir, trim_file, slurm_array_task_id, project, chemistry, self.nprocs, output_dir, os.path.basename(seq_dir), adapter_a, adapter_A, a_trim, h_filter, qiita_proj, final_output_dir))
        lines.append('echo "$cmd"')
        lines.append("date")
        lines.append("echo 'Executing: \"$cmd\"'")
        lines.append('eval "$cmd"')
        lines.append("date")

        file_path = os.path.join(seq_dir, 'Data/Fastq', project, file_base + '_qsub.sh')
        if os.path.exists(file_path):
            logging.error("%s is being overwritten." % file_path)

        with open(file_path, 'w') as f:
            for line in lines:
                f.write("%s\n" % line)

    def run(self, sample_sheet_path, seq_dir, slurm_array_task_id):
        metadata = self._get_run_config_info(sample_sheet_path)
        fastq_output = os.path.join(seq_dir, 'Data', 'Fastq')

        for read, project, forward_adapter, reverse_adapter, polyg_trimming, human_filtering, qiita_id in metadata:
            # extract the numerical component from the project string:
            # nnnnnn_ddddd
            proj_check = project.split('_')[-1]
            input_count = 0
            working_dir = os.path.join(fastq_output, project)
            for root, dirs, files in os.walk(working_dir):
                for some_file in files:
                    if some_file.endswith('.fastq.gz'):
                        if '_R1_' in some_file:
                            input_count += 1

            split_count = 1
            if input_count > 2000:
                split_count = 16
            elif input_count <= 2000 and input_count > 1000:
                split_count = 10
            elif input_count <= 1000 and input_count > 500:
                split_count = 4

            trim_file = 'split_file_'

            job_count = split_count + 1
            file_base = os.path.basename(self.seq_dir)

            for root, dirs, files in os.walk(working_dir):
                for some_file in files:
                    if trim_file in some_file:
                        some_path = os.path.join(root, some_file)
                        logging.debug("deleting %s..." % some_path)
                        os.remove(some_path)

            # we're only supposed to consider the contents of the
            # working_directory itself, not the sub-directories,
            # but for now assume that there aren't any sub-dirs.
            file_base_file_list = []
            for root, dirs, files in os.walk(working_dir):
                for some_file in files:
                    if some_file.endswith('.fastq.gz'):
                        if '_R1_' in some_file:
                            file_base_file_list.append(some_file)

            line_count = len(file_base_file_list)
            lines_per_split = (line_count + split_count - 1) / split_count

            self._split_list(file_base_file_list, lines_per_split, trim_file, 'extension', working_dir)

            # TODO: Modify create_qsub_file to take the member variables directly.

            self.create_qsub_file(self.seq_dir,
                                    project,
                                    file_base,
                                    self.home,
                                    self.email_list,
                                    trim_file,
                                    slurm_array_task_id,
                                    self.chemistry,
                                    self.output_dir,
                                    forward_adapter,
                                    reverse_adapter,
                                    polyg_trimming,
                                    human_filtering,
                                    qiita_id,
                                    self.final_output_dir)

            cmd = ['sbatch',
                   '--qos=seq_proc',
                   '--parsable',
                   '--array=0-%d' % split_count - 1,
                   os.path.join(seq_dir, 'Data', 'Fastq', project, file_base + '_qsub.sh')]

            # TODO We may or may not want to wait
            pbs_job_id = self.execute_sbatch_job_and_wait(cmd)

            pbs_job_array = "--dependency=afterok:$pbs_job_id"

            proj_array = project

            # TODO FastQCJobs should be created from the Pipeline, rather than
            #  within another job.
            ### submit fastqc for current project directory
            # fastqc_process "${seqdir}" "${output_dir}" "${fastq_output}" "${pbs_job_id}" "${pbs_job_array}" "${proj_array[@]}"






