from sequence_processing_pipeline.Job import Job
from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join, basename, dirname
from os import walk, remove
import logging


class QCJob(Job):
    def __init__(self, root_dir, sample_sheet_path, fpmmp_path, nprocs, mmi_db_path):
        super().__init__()
        self.root_dir = root_dir
        metadata = self._process_sample_sheet(sample_sheet_path)
        self.project_data = metadata['projects']
        self.trim_file = 'split_file_'
        self.fpmmp_path = fpmmp_path
        self.nprocs = nprocs
        self.chemistry = metadata['chemistry']
        # self.job_script_path = join(self.root_dir, 'human-filtering.sh')
        self.mmi_db_path = mmi_db_path
        self.stdout_log_path = 'localhost:' + join(self.root_dir, 'human-filtering.out.log')
        self.stderr_log_path = 'localhost:' + join(self.root_dir, 'human-filtering.err.log')

    def run(self, queue_name, node_count, nprocs, wall_time_limit):
        metadata = []

        for project in self.project_data:
            fastq_files = self._find_fastq_files(project['project_name'])
            split_count = self._generate_split_count(len(fastq_files))
            self._clear_trim_files()
            lines_per_split = (len(fastq_files) + split_count - 1) / split_count
            trim_files = self._generate_trim_files(fastq_files, lines_per_split)
            if len(trim_files) != split_count:
                logging.warning("The number of trim files does not equal the number of files we were supposed to have.")

            script_path = self._generate_job_script(queue_name,
                                                    node_count,
                                                    nprocs,
                                                    wall_time_limit,
                                                    split_count,
                                                    project['project_name'],
                                                    project['adapter_a'],
                                                    project['adapter_A'],
                                                    project['a_trim'],
                                                    project['h_filter'],
                                                    project['qiita_proj'])

            d = {'project_name': project['project_name']}
            # pbs_job_array="--dependency=afterok:$pbs_job_id"

            # TODO: replace w/job_info information, etc.
            fastq_output = None
            pbs_job_id = None
            pbs_job_array = None
            proj_array = None

            d['job_info'] = self.qsub(script_path, None, None)
            d['fastqc_params'] = {}
            d['fastqc_params']['seqdir'] = self.root_dir
            # for now, define output_dir in two locations until the location
            # is settled.
            d['fastqc_params']['output_dir'] = join(self.root_dir, 'human-filtering_output')
            d['fastqc_params']['fastq_output'] = fastq_output
            d['fastqc_params']['pbs_job_id'] = pbs_job_id
            d['fastqc_params']['pbs_job_array'] = pbs_job_array
            d['fastqc_params']['proj_array'] = proj_array
            metadata.append(d)

        return metadata

    def _generate_trim_files(self, fastq_files, split_count):
        def _chunk_list(some_list, n):
            # taken from https://bit.ly/38S9O8Z
            for i in range(0, len(some_list), n):
                yield some_list[i:i + n]

        # put these split files in the same location as the fastq_files for
        # the project. Assume all filepaths in fastq_files have the same
        # result for dirname().
        destination_path = dirname(fastq_files[0])

        new_files = []
        count = 0
        for chunk in _chunk_list(fastq_files, split_count):
            trim_file_name = '%s%d' % (self.trim_file, count)
            trim_file_path = join(destination_path, trim_file_name)
            with open(trim_file_path, 'w') as f:
                for line in chunk:
                    f.write("%s\n" % line)
            new_files.append(trim_file_path)
            count += 1

        return new_files

    def _clear_trim_files(self):
        # remove all files with a name beginning in self.trim_file.
        # assume cleaning the entire root_dir is overkill, but won't
        # hurt anything.
        for root, dirs, files in walk(self.root_dir):
            for some_file in files:
                if self.trim_file in some_file:
                    some_path = join(root, some_file)
                    remove(some_path)

    def _process_sample_sheet(self, sample_sheet_path):
        # Some column headers in some sample sheets appear to use
        # different names for the same metadata. This dict maps the names
        # commonly found to the internal name used by the code.
        # TODO: Not certain BarcodesAreRC == PolyGTrimming == a_trim. Confirm
        name_map = {'Sample_Project': 'project_name', 'Project': 'project_name', 'QiitaID': 'qiita_proj',
                    'BarcodesAreRC': 'a_trim', 'PolyGTrimming': 'a_trim', 'ForwardAdapter': 'adapter_a',
                    'ReverseAdapter': 'adapter_A', 'HumanFiltering': 'h_filter'}

        sheet = KLSampleSheet(sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % sample_sheet_path
            raise PipelineError(s)

        header = valid_sheet.Header
        chemistry = header['chemistry']

        bioinformatics = valid_sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.

        l = []
        for i in range(0, len(bioinformatics)):
            l.append({})

        for item in bioinformatics:
            my_series = bioinformatics[item]
            for index, value in my_series.items():
                key = name_map[item]
                # convert all manner of positive/negative implying strings
                # into a true boolean value.
                if key in ['a_trim', 'h_filter']:
                    value = value.strip().lower()
                    value = True if value in ['true', 'yes'] else False
                l[index][key] = value

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry, 'projects': l}

    def _find_fastq_files(self, project_name):
        search_path = join(self.root_dir, 'Data', 'Fastq', project_name)
        l = []
        for root, dirs, files in walk(search_path):
            for some_file in files:
                some_file = some_file.decode('UTF-8')
                if some_file.endswith('fastq.gz'):
                    if '_R1_' in some_file:
                        some_path = join(self.root_dir, some_file)
                        l.append(some_path)
        return l

    def _generate_split_count(self, count):
        if count > 2000:
            return 16
        elif count <= 2000 and count > 1000:
            return 10
        elif count <= 1000 and count > 500:
            return 4

        return 1

    def _generate_job_script(self, queue_name, node_count, nprocs, wall_time_limit, split_count, project_name,
                             adapter_a, adapter_A, a_trim, h_filter, qiita_proj):
        lines = []

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying params
        # w/environment variables is to do the same on QSUB but with -v instead of
        # --export. The syntax of the values are otherwise the same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        # lines.append("#PBS -N {}".format(project_name))
        lines.append("#PBS -N %s" % project_name)

        # what torque calls a queue, slurm calls a partition
        # SBATCH -p SOMETHING -> PBS -q SOMETHING
        # (Torque doesn't appear to have a quality of service (-q) option. so this
        # will go unused in translation.)
        lines.append("#PBS -q %s" % queue_name)

        # request one node
        # Slurm --ntasks-per-node=<count> -> -l ppn=<count>	in Torque
        lines.append("#PBS -l nodes=%d:ppn=%d" % (node_count, nprocs))

        # Slurm --export=ALL -> Torque's -V
        lines.append("#PBS -V")

        # Slurm walltime limit --time=24:00:00 -> Torque's -l walltime=<hh:mm:ss>
        # using the larger value found in the two scripts (72 vs ? hours)
        lines.append("#PBS -l walltime=%d:00:00" % wall_time_limit)

        # send email to charlie when a job starts and when it terminates or
        # aborts. This is used to confirm the package's own reporting
        # mechanism is reporting correctly.
        lines.append("#PBS -m bea")

        # specify your email address
        lines.append("#PBS -M ccowart@ucsd.edu")

        # min mem per CPU: --mem-per-cpu=<memory> -> -l pmem=<limit>
        # taking the larger of both values (10G > 6G)
        # lines.append("#PBS -l pmem=10gb")
        # revisit the mem stuff later.

        # --output -> -o
        lines.append("#PBS -o %s" % self.stdout_log_path)
        lines.append("#PBS -e %s" % self.stderr_log_path)

        # array input files are labeled file0, file1,...filen-1
        # TODO: watch edge-case where only one split_count file
        #  is generated. This line will become '-t 0-0'. I believe
        #  that is valid, but it may be better to state as '-t 0'.
        lines.append("#PBS -t 0-%d" % split_count - 1)

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1

        # We probably do not need to activate this Python environment, but I
        # will store it here in comments.
        # source ~/miniconda3/bin/activate test_env_2

        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use root_dir instead.
        lines.append("set -x")
        lines.append("cd %s" % self.root_dir)
        lines.append("date '+%s' > human-filtering.job.log")
        lines.append("module load fastp_0.20.1 samtools_1.12 minimap2_2.18")

        # for now, set the output from fpmmp.sh to be in two separate distinct
        # directories until we can confirm its usual location. We will make
        # the directories to ensure they are fresh.
        lines.append("rm -rf human-filtering_output")
        lines.append("mkdir human-filtering_output")
        lines.append("rm -rf human-filtering_final_output")
        lines.append("mkdir human-filtering_final_output")
        output_dir = join(self.root_dir, 'human-filtering_output')
        final_output_dir = join(self.root_dir, 'human-filtering_final_output')

        cmd = []
        cmd.append(self.fpmmp_path)
        cmd.append('-d %s' % self.root_dir)
        cmd.append('-D %s/Data/Fastq' % self.root_dir)
        # even though we have a list of all of the trim_files0-(n-1),
        # we use this syntax so that we only need to qsub one job file, insteaf
        # of n job files.
        cmd.append('-S %s${PBS_ARRAYID}' % self.trim_file)
        # cmd.append('-p %s' % project_name)
        cmd.append('-p %s' % project_name)
        cmd.append('-x %s' % self.mmi_db_path)
        # cmd.append('-C %s' % chemistry)
        cmd.append('-C %s' % self.chemistry)
        cmd.append('-c %s' % self.nprocs)
        cmd.append('-o %s' % output_dir)
        cmd.append('-O %s' % basename(self.root_dir))
        cmd.append('-a %s' % adapter_a)
        cmd.append('-A %s' % adapter_A)
        cmd.append('-g %s' % a_trim)
        cmd.append('-G %s' % h_filter)
        cmd.append('-q %s' % qiita_proj)
        cmd.append('-f %s' % final_output_dir)
        lines.append(' '.join(cmd))

        # unlike w/BCL2FASTQJob, multiple human-filtering.sh scripts
        # will be generated, one for each project defined in the
        # sample sheet.
        job_script_path = join(self.root_dir, 'human-filtering-%s.sh' % project_name)
        with open(job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                f.write("%s\n" % line)


if __name__ == '__main__':
    import json

    hf2job = HumanFilterJob('.', './good-sample-sheet.csv', '.', './fpmmp.sh', 16)
    results = hf2job._process_sample_sheet('./good-sample-sheet.csv')
    for key in results:
        print(key)
        print(json.dumps(results[key], indent=2))
