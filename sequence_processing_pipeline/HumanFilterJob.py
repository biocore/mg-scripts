from os import makedirs, listdir
from os.path import join, dirname, abspath, exists, isfile, basename
from sequence_processing_pipeline.TorqueJob import TorqueJob
from time import sleep
from zipfile import ZipFile
import logging
import os
import re
import shutil


#source ~/.seqvars
#module load fastp_0.20.1 samtools_1.12 minimap2_2.18


class HumanFilterJob(TorqueJob):
    def __init__(self, root_dir, project, sample_sheet_path, output_directory, output2_directory, forward, reverse, project2=None):
        super().__init__()
        logging.debug("HumanFilterJob Constructor called")

        self.job_script_path = join(self.project_fastq_dir, 'fastqc_qsub.sh')
        self.job_name = "%s_%s_%s" % (project, forward, reverse)
        self.nprocs = 16
        self.sample_sheet_path = sample_sheet_path
        self.stdout_log_path = join(self.project_fastq_dir, 'something_%s_%s.out.log' % forward, reverse)
        self.stderr_log_path = join(self.project_fastq_dir, 'something_%s_%s.err.log' % forward, reverse)
        self.output_dir_path = abspath(output_directory)
        # sometimes called final_output_directory and used by CMI.
        self.output2_dir_path = abspath(output2_directory)
        self.bcl2fastq_bin_path = "~/bcl2fastq_2_20/bin/bcl2fastq"

        # this should be written out to
        # $final_output/${project}/empty_file_list.txt if needed.
        self.empty_file_list_txt = []

        self.trim_list = self._generate_trim_list()
        self._create_directories(root_dir, project, project2)

    def _generate_trim_list(self):
        l = []
        for item in listdir(self.project_fastq_dir):
            if isfile(join(self.project_fastq_dir, item)):
                if item.endswith('.fastq.gz'):
                    if '_R1_' in item:
                        l.append(item)
        return l

    def _create_directories(self, root_dir, project, project2=None):
        # this Job creates a lot of new directories. It's easier to simply
        # create all of the ones needed and walk the tree after processing
        # to remove empty directories. What's important for the copy and move
        # operations is that the destination directories exist.
        path_list = []

        self.project = project
        self.root_dir = abspath(root_dir)
        self.project_dir = join(self.root_dir, project)
        self.fastq_dir = join(self.root_dir, 'Data', 'Fastq')
        self.project_fastq_dir = join(self.fastq_dir, project)

        path_list.append(self.root_dir)
        path_list.append(self.project_dir)
        path_list.append(self.fastq_dir)
        path_list.append(self.project_fastq_dir)

        path_list.append(join(self.project_fastq_dir, 'html'))
        path_list.append(join(self.project_fastq_dir, 'json'))

        path_list.append(join(self.project_dir, 'fastp_qc', 'fastp_fastqc'))
        path_list.append(join(self.project_dir, 'fastp_qc', 'fastp_logs'))

        filtered_sequences_location_1 = join(self.project_dir, 'filtered_sequences')
        path_list.append(filtered_sequences_location_1)
        path_list.append(join(filtered_sequences_location_1, 'trim_logs'))
        path_list.append(join(filtered_sequences_location_1, 'zero_files'))

        trimmed_sequences_location_1 = join(self.project_dir, 'trimmed_sequences')
        path_list.append(trimmed_sequences_location_1)

        if project2:
            filtered_sequences_location_2 = join(self.project_dir, project2, 'filtered_sequences')
            path_list.append(filtered_sequences_location_2)
            path_list.append(join(filtered_sequences_location_2, 'trim_logs'))
            path_list.append(join(filtered_sequences_location_2, 'zero_files'))

            trimmed_sequences_location_2 = join(self.project_dir, project2, 'trimmed_sequences')
            path_list.append(trimmed_sequences_location_2)

        path_list.append(join(self.root_dir, 'html'))
        path_list.append(join(self.project_dir, 'html'))

        path_list.append(join(self.root_dir, 'json'))
        path_list.append(join(self.project_dir, 'json'))

        path_list.append(join(self.root_dir, 'trim_logs'))
        path_list.append(join(self.project_dir, 'trim_logs'))

        path_list.append(join(self.root_dir, 'zero_files'))
        path_list.append(join(self.project_dir, 'zero_files'))

        path_list.append(join(self.project_dir, 'fastp_qc', 'fastp_fastqc'))
        path_list.append(join(self.project_dir, 'fastp_qc', 'fastp_logs'))
        path_list.append(join(self.project_dir, 'fastp_qc', 'processed_qc'))

        path_list.append(join(self.project_fastq_dir, 'fastp_qc', 'fastp_fastqc'))
        path_list.append(join(self.project_fastq_dir, 'fastp_qc', 'fastp_logs'))
        path_list.append(join(self.project_fastq_dir, 'fastp_qc', 'processed_qc'))

        path_list.append(join(self.output2_dir_path, self.project, 'amplicon'))

        # intentionally wrong location.
        path_list.append(join(self.project_dir, 'wrong_location'))

        for path in path_list:
            makedirs(path, exist_ok=True)

    def zip_files(self, file_paths, zip_path):
        with ZipFile(zip_path, 'w') as zip:
            for file in file_paths:
                zip.write(file)

    def copy_directory_contents(self, source_directory, destination_directory, ends_with=None, contains=None,
                                starts_with=None, delete_after=False, uid=None, gid=None):
        if uid and not gid:
            raise ValueError("You must define both uid and gid")
        elif not gid and uid:
            raise ValueError("You must define both uid and gid")

        files = listdir(source_directory)
        # convert from bytes output of listdir to string
        files = [str(x) for x in files]
        for some_file in files:
            match_count = 0
            count = 0
            if starts_with:
                count += 1
                if some_file.startswith(starts_with):
                    match_count += 1
            elif ends_with:
                count += 1
                if some_file.endswith(ends_with):
                    match_count += 1
            elif contains:
                count += 1
                if contains in some_file:
                    match_count += 1

            if count == match_count and count != 0:
                # if multiple filters were provided, we need to match all of
                # them before moving/copying files.
                if delete_after:
                    shutil.move(join(source_directory, some_file), destination_directory)
                else:
                    shutil.copyfile(join(source_directory, some_file), destination_directory)

                if uid and gid:
                    os.chown(join(destination_directory, some_file), uid, gid)

    def possible_amplicon(self, adapter_a, adapter_A, chemistry, qiita_proj, copy_file):
        for file in self.trim_list:
            filepath1 = file
            filepath2 = file.replace('_R1_00', '_R2_00')
            filename1 = basename(file)
            filename2 = filename1.replace('_R1_00', '_R2_00')

            filename1_short = filename1.replace('.fastq.gz', '')
            filename2_short = filename2.replace('.fastq.gz', '')

            filename_index1 = filename1.replace('_R1_', '_I1_')
            filename_index2 = filename2.replace('_R1_', '_I1_')

            cmd_list = ['fastp']

            if adapter_a != None and adapter_A != None:
                cmd_list.append('--adapter_sequence')
                cmd_list.append(adapter_a)
                cmd_list.append('--adapter_sequence_r2')
                cmd_list.append(adapter_A)

            f1 = join(self.project_dir, 'trimmed_sequences', filename1_short + '.fastp.fastq.gz')
            f2 = join(self.project_dir, 'trimmed_sequences', filename2_short + '.fastp.fastq.gz')

            cmd_list += ['-l', '100',
                         '-i', filepath1,
                         '-I', filepath2,
                         '-w', self.nprocs,
                         '-j', join(self.project_dir, 'json', filename1_short + '.json'),
                         '-h', join(self.project_dir, 'html', filename1_short + '.html'),
                         '-o', f1,
                         '-O', f2,
                         '-R', filename1_short + '_report'
                         ]

            out, err, rc = self._system_call(cmd_list)

            size1 = os.stat(f1).st_size
            size2 = os.stat(f2).st_size

            ### amplicon/16s data. include index files
            if chemistry == 'Amplicon':
                dst = join(self.project_dir, 'trimmed_sequences')
                self.copy_directory_contents(self.project_dir, dst, contains=filename_index1)
                self.copy_directory_contents(self.project_dir, dst, contains=filename_index2)

            if size1 <= 500 or size2 <= 500:
                src = join(self.project_dir, 'trimmed_sequences')
                dst = join(self.project_dir, 'zero_files')
                self.copy_directory_contents(src, dst, contains=filename1_short, delete_after=True)
                self.copy_directory_contents(src, dst, contains=filename2_short, delete_after=True)

                a = join(src, filename1_short + '.fastp.fastq.gz')
                b = join(src, filename2_short + '.fastp.fastq.gz')

                self.empty_file_list_txt.append(a)
                self.empty_file_list_txt.append(size1)
                self.empty_file_list_txt.append(b)
                self.empty_file_list_txt.append(size2)

        for some_file in self.trim_list:
            transfer1 = some_file.replace('.fastq.gz', '')
            transfer2 = transfer1.replace('_R1_00', '_R2_00')

            if chemistry == 'Amplicon':
                src = join(self.project_dir, 'trimmed_sequences')
                # unsure what dst should be.
                dst = join(self.project_dir, 'wrong_location')
                self.copy_directory_contents(src, dst, ends_with='.fastp.fastq.gz')
            else:
                # unsure what dst should be.
                dst = join(self.project_dir, 'wrong_location')
                self.copy_directory_contents(self.project_dir, dst, ends_with=some_file)
                filter = transfer2 + '.fastq.gz'
                self.copy_directory_contents(self.project_dir, dst, ends_with=filter)

            if qiita_proj == 'NA':
                if chemistry == 'Amplicon' and copy_file == 'TRUE':
                    src = join(self.project_dir, 'trimmed_sequences')
                    dst = '/qmounts/qiita_data/uploads/%s/' % qiita_proj
                    self.copy_directory_contents(src, dst, ends_with='.fastq.gz', contains='_R1_', uid='5500', gid='5500')
                    self.copy_directory_contents(src, dst, ends_with='.fastq.gz', contains='_R2_', uid='5500', gid='5500')
                    self.copy_directory_contents(src, dst, ends_with='.fastq.gz', contains='_I1_', uid='5500', gid='5500')
                elif copy_file == 'TRUE':
                    dst = '/qmounts/qiita_data/uploads/%s/' % qiita_proj
                    src1 = join(self.project_dir, 'trimmed_sequences', transfer1 + '.fastp.fastq.gz')
                    self.copy_directory_contents(src1, dst, uid='5500', gid='5500')
                    src2 = join(self.project_dir, 'trimmed_sequences', transfer2 + '.fastp.fastq.gz')
                    self.copy_directory_contents(src2, dst, uid='5500', gid='5500')

        # commented-out: Not  sure why we would have to do this
        # rsync -avp --progress ${final_output}/${project}/ ${final_output_dir}/${base_seq_dir}

    def possible_amplicon2(self, adapter_a, adapter_A, final_output, copy_file, qiita_proj):
        partial_path = join(final_output, self.project)
        for file in self.trim_list:
            parent_dir = dirname(file)
            tmp = parent_dir.split('/')
            # if tmp begins with a leading empty string '', this will remove it.
            tmp = [x for x in tmp if x]
            project_dir = tmp[0]
            filename1 = basename(file)
            filename1_short = filename1.replace('.fastq.gz', '')
            filename2 = filename1.replace('_R1_00', '_R2_00')
            filename2_short = filename2.replace('.fastq.gz', '')
            if adapter_a == 'NA' or adapter_A == 'NA':
                pass
            # $fastp
            # -l
            # 100
            # -i
            # $filename1
            # -I
            # $filename2
            # -w
            # $NPROCS
            # --stdout
            # -j
            # ${dir}/${project}/json/${filename1_short}.json
            # -h
            # ${dir}/${project}/html/${filename1_short}.html
            # | $minimap2 -ax sr -t $NPROCS $db - -a
            # | $samtools fastq -@ $NPROCS -f 12 -F 256 -1 $final_output/${project}/filtered_sequences/${filename1_short}.trimmed.fastq.gz -2 $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz
            else:
                pass
            # $fastp
            # --adapter_sequence
            # ${adapter_a}
            # --adapter_sequence_r2
            # ${adapter_A}
            # -l
            # 100
            # -i $filename1 -I $filename2 -w $NPROCS --stdout -j ${dir}/${project}/json/${filename1_short}.json -h ${dir}/${project}/html/${filename1_short}.html
            # | $minimap2 -ax sr -t $NPROCS $db - -a
            # | $samtools fastq -@ $NPROCS -f 12 -F 256 -1 $final_output/${project}/filtered_sequences/${filename1_short}.trimmed.fastq.gz -2 $final_output/${project}/filtered_sequences/${filename2_short}.trimmed.fastq.gz

###
            cmd_list = ['fastp']

            if adapter_a != None and adapter_A != None:
                cmd_list.append('--adapter_sequence')
                cmd_list.append(adapter_a)
                cmd_list.append('--adapter_sequence_r2')
                cmd_list.append(adapter_A)

            f1 = join(self.project_dir, 'trimmed_sequences', filename1_short + '.fastp.fastq.gz')
            f2 = join(self.project_dir, 'trimmed_sequences', filename2_short + '.fastp.fastq.gz')

            cmd_list += ['-l', '100',
                         '-i', filepath1,
                         '-I', filepath2,
                         '-w', self.nprocs,
                         '-j', join(self.project_dir, 'json', filename1_short + '.json'),
                         '-h', join(self.project_dir, 'html', filename1_short + '.html'),
                         '-o', f1,
                         '-O', f2,
                         '-R', filename1_short + '_report'
                         ]
###
            f1 = join(partial_path, 'filtered_sequences', filename1_short + '.trimmed.fastq.gz')
            f2 = join(partial_path, 'filtered_sequences', filename2_short + '.trimmed.fastq.gz')

            size1 = os.stat(f1).st_size
            size2 = os.stat(f2).st_size

            if size1 <= 500 or size2 <= 500:
                src = join(partial_path, 'trimmed_sequences')
                dst = join(partial_path, 'zero_files')

                # mv $final_output/${project}/filtered_sequences/${filename1_short}* $final_output/${project}/filtered_sequences/${filename2_short}* ${final_output}/${project}/zero_files
                self.copy_directory_contents(self.dir, dst, contains=filename1_short, delete_after=True)

                self.empty_file_list_txt.append(
                    join(partial_path, 'filtered_sequences$', filename1_short + '.trimmed.fastq.gz'))
                self.empty_file_list_txt.append(size1)
                self.empty_file_list_txt.append(
                    join(partial_path, 'filtered_sequences$', filename2_short + '.trimmed.fastq.gz'))
                self.empty_file_list_txt.append(size2)

        ### always executes
        if copy_file:
            for file_transfer in self.trim_list:
                # transfer_1=$(echo "$file_transfer" | cut -f1 -d".")
                # transfer_2=$(echo "$transfer_1" | sed -e 's/_R1_00/_R2_00/g')
                # rsync -avp --progress ${final_output}/${project}/filtered_sequences/${transfer_1}.trimmed.fastq.gz ${final_output}/${project}/filtered_sequences/${transfer_2}.trimmed.fastq.gz ${final_output_dir}/${base_seq_dir}/${project}/filtered_sequences/
                if qiita_proj != 'NA' and copy_file == 'TRUE':
                    pass

        # sudo -u qiita rsync -avp --chown 5500:5500 ${final_output}/${project}/filtered_sequences/${transfer_1}.trimmed.fastq.gz ${final_output}/${project}/filtered_sequences/${transfer_2}.trimmed.fastq.gz /qmounts/qiita_data/uploads/${qiita_proj}/

        # rsync -avp --progress ${final_output} ${final_output_dir}/

    def run(self, rundir, filter_db, a_trim, h_filter, chemistry, qiita_proj, copy_file, adapter_a, adapter_A, nprocs,
             final_output):
        '''
        # we may actually need to make the job script in a sub-dir.
        self._make_job_script()
        job_info = self.qsub(self.job_script_path, None, None)

        # if the job returned successfully, we may want to send an
        # email to the user. If an error occurs, the user will
        # be notified by the email triggered by PipelineError().
        '''
        self._make_job_script()

        if not isinstance(a_trim, bool):
            raise ValueError("a_trim must be a boolean value.")

        if not isinstance(h_filter, bool):
            raise ValueError("h_filter must be a boolean value.")

        if exists(join(rundir, 'run_config.h')):
            pass
        # TODO: source the file - we'll need to find an alternative method
        #  to take in the settings

        filter_db = "/databases/bowtie/Human/Human"

        # move sample sheets into place
        # pull everything from original location to top level output
        self.copy_directory_contents(self.dir, self.output_dir, '.csv')

        ### amplicon MiSeq sequencing
        if not a_trim and not h_filter:
            ### just copy data to qiita study.  16s data.
            if chemistry == 'Amplicon':
                src = join(self.dir, self.project)
                dst = join(self.output2_dir_path, self.project, 'amplicon')
                self.copy_directory_contents(src, dst, '.fastq.gz')
            if qiita_proj:
                if qiita_proj != 'NA' and copy_file != 'FALSE':
                    src = join(self.dir, self.project)
                    dst = join('qmounts', 'qiita_data', 'uploads', qiita_proj)
                    # these probably need to change ownership to qiita and group to the proper group as well.
                    self.copy_directory_contents(src, dst, '.fastq.gz')

        if a_trim and not h_filter:
            self.possible_amplicon(adapter_a, adapter_A, chemistry, qiita_proj, copy_file)
        elif a_trim and h_filter:
            self.possible_amplicon2(adapter_a, adapter_A, final_output, copy_file, qiita_proj)

        if not h_filter:
            # if the files failed to transfer then exit with -1
            pass
        else:
            exit(0)

        do_bowtie = True

        if do_bowtie:
            for some_file in self.trim_list:
                final_output = join(self.dir, self.project, 'filtered_sequences')
                parent_dir = dirname(some_file)
                tmp = parent_dir.split('/')
                # if tmp begins with a leading empty string '', this will remove it.
                tmp = [x for x in tmp if x]
                project_dir = tmp[0]
                filename1 = basename(some_file)
                filename1_short = filename1 + '.fastq.gz'
                filename2 = filename1.replace('_R1_', '_R2_')
                filename2_short = filename2 + '.fastq.gz'

                from subprocess import Popen, PIPE

                p1 = join(self.fastp_qc_output, filename1)
                p2 = join(self.fastp_qc_output, filename2)

                bowtie_log = join(final_output, 'trim_logs', filename1_short + '.log')
                with open(bowtie_log, 'a') as f:
                    bowtie = Popen(['bowtime', '-p', nprocs, '-x', filter_db, '--very-sensitive', '-1', p1, '-2', p2],
                                   stdout=PIPE, stderr=f)
                    samtools1 = Popen(['samtools', 'view', '-f', '12', '-F', '256'], stdout=PIPE, stderr=PIPE,
                                      stdin=bowtie.stdout)
                    samtools2 = Popen(['samtools', 'sort', '-@', '16', '-n'], stdout=PIPE, stderr=PIPE,
                                      stdin=samtools1.stdout)
                    samtools3 = Popen(['samtools', 'view', '-bS'], stdout=PIPE, stderr=PIPE, stdin=samtools2.stdout)
                    bedtools = Popen(['bedtools', 'bamtofastq', '-i', '-', '-fq',
                                      join(final_output, filename1_short + '.trimmed.fastq'), '-fq2',
                                      join(final_output, filename1_short + '.trimmed.fastq')], stdout=f, stderr=PIPE,
                                     stdin=samtools3.stdout)

                    stdout, stderr = bedtools.communicate()
                    return_code = bedtools.returncode

                    file1 = join(final_output, filename1_short + '.trimmed.fastq')
                    self.zip_files(file1, file1 + '.zip')
                    file2 = join(final_output, filename2_short + '.trimmed.fastq')
                    self.zip_files(file2, file2 + '.zip')

                    # TODO: Find an alternate method to check and see if the two
                    #  files gzipped above should be added to zero_files.txt and
                    #  the gzips deleted.
                    for some_file in self.trim_list:
                        tmp = some_file.split('_')
                        # extract 'X00178928_S1_L003' from
                        # 'X00178928_S1_L003_R1_001.fastq.gz' ...
                        tmp = '_'.join(tmp[0:3])

                        # this code may be vestigial. the contents of both
                        # conditionals were commented out.
                        if qiita_proj != 'NA':
                            # echo copying files ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz
                            pass
                        else:
                            # echo syncing file ${final_output}/${file}_R1_*.trimmed.fastq.gz ${final_output}/${file}_R2_*.trimmed.fastq.gz
                            pass

    def _make_job_script(self):
        lines = []

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying params
        # w/environment variables is to do the same on QSUB but with -v instead of
        # --export. The syntax of the values are otherwise the same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        lines.append("#PBS -N %s" % self.job_name)

        # what torque calls a queue, slurm calls a partition
        # SBATCH -p SOMETHING -> PBS -q SOMETHING
        # (Torque doesn't appear to have a quality of service (-q) option. so this
        # will go unused in translation.)
        lines.append("#PBS -q long8gb")

        # request one node
        # Slurm --ntasks-per-node=<count> -> -l ppn=<count>	in Torque
        lines.append("#PBS -l nodes=1:ppn=%d" % self.nprocs)

        # Slurm --export=ALL -> Torque's -V
        lines.append("#PBS -V")

        # Slurm walltime limit --time=24:00:00 -> Torque's -l walltime=<hh:mm:ss>
        # using the larger value found in the two scripts (72 vs ? hours)
        lines.append("#PBS -l walltime=72:00:00")

        # send email to ccowart when a job starts and when it terminates or
        # aborts. This is used to confirm the package's own reporting
        # mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # specify your email address
        lines.append("#PBS -M ccowart@ucsd.edu")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        #lines.append("#PBS -l pmem=10gb")
        # revisit the mem stuff later.

        # --output -> -o
        lines.append("#PBS -o %s" % self.stdout_log_path)
        lines.append("#PBS -e %s" % self.stderr_log_path)

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1

        # We probably do not need to activate this Python environment, but I
        # will store it here in comments.
        # source ~/miniconda3/bin/activate test_env_2

        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use root_dir instead.
        lines.append("cd %s"  % self.root_dir)
        '''
                lines.append(
                    'cmd="sh ~/seq_proc_dev/fpmmp.sh -d %s -D /Data/Fastq -S %s\%s -p %s -C %s -c %d -o %s -O %s -a %s -A %s -g %s -G %s -q %s -f %s"' % (
                    seq_dir, trim_file, slurm_array_task_id, project, chemistry, self.nprocs, output_dir,
                    os.path.basename(seq_dir), adapter_a, adapter_A, a_trim, h_filter, qiita_proj, final_output_dir))
        '''
        lines.append('cmd="%s --sample-sheet %s --mask-short-adapter-reads \
                      1 -R . -o %s --loading-threads 8 --processing-threads \
                      8 --writing-threads 2 --create-fastq-for-index-reads \
                      --ignore-missing-bcls"' % (self.bcl2fastq_bin_path,
                                                 self.sample_sheet_path,
                                                 self.output_dir_path))
        lines.append("echo $cmd")
        lines.append("eval $cmd")
        lines.append("return_code=$?")
        p = "%s/%s.return_code" % (self.root_dir, self.job_name)
        lines.append("echo $returncode > %s" % p)
        lines.append("date >> %s" % p)

        with open(self.job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % self.job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                line = re.sub('\s+', ' ', line)
                f.write("%s\n" % line)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    job = HumanFilterJob('./good-bcl-directory',
                       'THDMI_US_99999',
                       './good-bcl-directory/good-sample-sheet.csv',
                       './good-bcl-directory/Data/Fastq/Output')
    job.run()


