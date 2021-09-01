import os
import logging
from exceptions import PipelineError
import re


class HumanFilter:
    def __init__(self, SequenceDirectoryObject, nprocs):
        self.sdo = SequenceDirectoryObject
        self.seq_dir = self.sdo.seq_dir
        self.csv_files = []
        self.nprocs = nprocs

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

    def _generate_run_config_file(self):
        skip_these_files = ['sav.csv']
        results = []

        for csv_file in self.csv_files:
            if csv_file in skip_these_files:
                logging.debug("Skipping %s..." % csv_file)
            else:
                with open(csv_file, 'r') as f:
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
                    # the original bash scripts.)
                    l = list(set(metadata))
                    l.sort()

                    # write the sanitized data out to legacy file.
                    # with open(self.run_config_file_path, 'w') as f2:

                    # write the sanitized data out to a file that can be
                    # associated with each csv_file.
                    run_config_file_path = '%s.run_config' % csv_file
                    with open(run_config_file_path, 'w') as f2:
                        for line in metadata:
                            # end lines w/proper UNIX-style newline, unless
                            # Windows '\r\n' is preferred.
                            f2.write("%s\n" % line)
                        results.append(run_config_file_path)
        if results:
            return results

    def qsub_file_generator(self, seq_dir, project, file_base, home_dir,
                            email_list, trim_file, slurm_array_task_id,
                            chemistry, output_dir, adapter_a, adapter_A,
                            a_trim, h_filter, qiita_proj, final_output_dir,
                            cmd):
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
        lines.append("#SBATCH --output={}}/filter_jobs/%x_%A_%a.out".format(home_dir))
        lines.append("#SBATCH --error={}}/filter_jobs/%x_%A_%a.err".format(home_dir))
        lines.append("#SBATCH --mail-type=ALL")
        lines.append("#SBATCH --mail-user=%s" % ','.join(email_list))
        # assume for now this hardcoded path is acceptable
        lines.append("source ~/miniconda3/bin/activate test_env_2")
        # assume for now the backslash below is intentional
        lines.append("file=%s\%s" % (trim_file, slurm_array_task_id))
        lines.append(
            "cmd=\"sh ~/seq_proc_dev/fpmmp.sh -d %s -D /Data/Fastq -S %s\%s -p %s -C %s -c %d -o %s -O %s -a %s -A %s -g %s -G %s -q %s -f %s\"" % (
                seq_dir, trim_file, slurm_array_task_id, project, chemistry, self.nprocs, output_dir, os.path.basename(seq_dir), adapter_a, adapter_A, a_trim, h_filter, qiita_proj, final_output_dir))
        lines.append("echo \"%s\"" % cmd)
        lines.append("date")
        lines.append("echo 'Executing: %s'" % cmd)
        lines.append("eval \"%s\"" % cmd)
        lines.append("date")

        file_path = os.path.join(seq_dir, 'Data/Fastq', project, file_base + '_qsub.sh')
        if os.path.exists(file_path):
            logging.error("%s is being overwritten." %  file_path)

        with open(file_path, 'w') as f:
            for line in lines:
                f.write("%s\n" % line)

    def _get_split_count(self, some_path):
        input_count = 0
        for root, dirs, files in os.walk(some_path):
            for some_file in files:
                results = re.search('.*R1.*\.fastq\.gz', some_file)
                if results:
                    input_count += 1

        if input_count > 2000:
            split_count = 16
        elif input_count <= 2000 and input_count > 1000:
            split_count = 10
        elif input_count <= 1000 and input_count > 500:
            split_count = 4
        else:
            split_count = 1

        return split_count

    def _process_single_run_config_file(self, run_config_file_path):
        with open(run_config_file_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                try:
                    project,\
                    forward_adapter,\
                    reverse_adapter,\
                    polyg_trimming,\
                    human_filtering,\
                    qiita_id = line.split(',')
                    fastq_output = os.path.join(self.seq_dir, 'Data', 'Fastq')
                    split_count = self.get_split_count(os.path.join(fastq_output, project))
                    trim_file = 'split_file_'



                    job_count = split_count + 1
                    file_base = os.path.basename(self.seq_dir)

                    # delete all files with 'split_file_' in them

                    # I believe this is the right subdirectory
                    mylist = []
                    for root, dirs, files in os.walk(os.path.join(fastq_output, project)):
                        for some_file in files:
                            if some_file.endswith('.fastq.gz'):
                                if '_R1_' in some_file:
                                    mylist.append(some_file)

                    # since we only use <filebase>_file_list.txt to split into n files,
                    # we don't need to actually create it.
                    line_count = len(mylist)
                    lines_per_split = (line_count + split_count - 1) / split_count

                    #  I believe fastq_output is the correct directory to output these files.
                    split_files = self._split_file(mylist, lines_per_split, 'split_file_', '', fastq_output)

                    self.qsub_file_generator(self.seq_dir,
                                             project,
                                             file_base,
                                             self.nprocs,
                                             home,
                                             email_list,
                                             trim_file,
                                             slurm_array_task_id,
                                             chemistry,
                                             output_dir,
                                             base_name,
                                             forward_adapter,
                                             reverse_adapter,
                                             polyg_trimming,
                                             human_filtering,
                                             qiita_id,
                                             final_output_dir,
                                             cmd)

                    # submit fastqc process  for current  directory.
                    # fastqc_process "${seqdir}" "${output_dir}" "${fastq_output}" "${pbs_job_id}" "${pbs_job_array}" "${proj_array[@]}"

                except ValueError as e:
                    raise PipelineError("%s does not contain the proper number of columns." % self.run_config_file_path)


    def run(self):
        # if not os.path.exists(self.run_config_file_path):
        # # run_param doesn't run anything, it actually generates the run_config.txt file. :)
        # self._run_param()

        # for now, let's assume this method needs to be run each time. There
        # doesn't appear to be a negative for re-doing the data. We are also
        # generating multiple files now - one for each csv file found. Also,
        # the name of the method has changed from _run_param() to
        # _generate_run_config_file().
        run_config_paths = self._generate_run_config_file()

        for run_config_path in run_config_paths:
            self._process_single_run_config_file(run_config_path)
