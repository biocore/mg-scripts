from sequence_processing_pipeline.TorqueJob import TorqueJob
from metapool import KLSampleSheet, validate_sample_sheet
from sequence_processing_pipeline.PipelineError import PipelineError
from os.path import join, basename, dirname
from os import walk, remove
import logging


class HumanFilter2Job(TorqueJob):
    def __init__(self, root_dir, sample_sheet_path, log_directory, fpmmp_path, nprocs):
        super().__init__()
        self.root_dir = root_dir
        metadata = self._process_sample_sheet(sample_sheet_path)
        self.project_data = metadata['projects']
        self.trim_file = 'split_file_'
        self.log_directory = log_directory
        self.fpmmp_path = fpmmp_path
        self.nprocs = nprocs
        self.chemistry = metadata['chemistry']

    def run(self):
        for project in self.project_data:
            fastq_files = self._find_fastq_files(project['project_name'])
            split_count = self._generate_split_count(len(fastq_files))
            job_count = split_count + 1
            file_base = basename(self.root_dir)  # this might need to be data/fastq instead
            self._clear_trim_files()
            self.lines_per_split = (len(fastq_files) + split_count - 1) / split_count
            trim_files = self._generate_trim_files(fastq_files, split_count)
            proj_array = []

            # TODO: Replace with appropriate relative paths
            output_dir = 'MyOutputDirectory'
            final_output_dir = 'MyFinalOutputDirectory'

            script_path = self._make_job_script(project['project_name'],
                                                self.chemistry,
                                                output_dir,
                                                project['adapter_a'],
                                                project['adapter_A'],
                                                project['a_trim'],
                                                project['h_filter'],
                                                project['qiita_proj'],
                                                final_output_dir)

            # TODO: We need to supply qsub with the Torque version of
            #  pbs_job_id=$(sbatch --qos=seq_proc --parsable --array=0-$(($split_count - 1)) ${seqdir}/Data/Fastq/${project}/${file_base}_qsub.sh)
            #  pbs_job_array="--dependency=afterok:$pbs_job_id"
            job_info = self.qsub(script_path, None, None)

            proj_array.append(project['project_name'])
            # fastqc_process "${seqdir}" "${output_dir}" "${fastq_output}" "${pbs_job_id}" "${pbs_job_array}" "${proj_array[@]}"

        # Don't bother creating uniq_run_config.txt at this time.
        # done < ${seqdir} / uniq_run_config.txt

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
                    f.write(line)
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
        valid_sheet = validate_sample_sheet(sheet)

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

    def _make_job_script(self, project_name, chemistry, output_dir, adapter_a, adapter_A, a_trim, h_filter, qiita_proj,
                         final_output_dir):
        lines = []

        lines.append("#!/bin/bash")
        # The Torque equiv for calling SBATCH on this script and supplying params
        # w/environment variables is to do the same on QSUB but with -v instead of
        # --export. The syntax of the values are otherwise the same.
        # -v <variable[=value][,variable2=value2[,...]]>

        # declare a name for this job to be sample_job
        lines.append("#PBS -N {}_%A_%a".format(project_name))

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
        lines.append("#PBS -o Localhost:{}/filter_jobs/%x_%A_%a.out".format(self.log_directory))
        lines.append("#PBS -e Localhost:{}/filter_jobs/%x_%A_%a.err".format(self.log_directory))

        # there is no equivalent for this in Torque, I believe
        # --cpus-per-task 1

        # We probably do not need to activate this Python environment, but I
        # will store it here in comments.
        # source ~/miniconda3/bin/activate test_env_2

        # By default, PBS scripts execute in your home directory, not the
        # directory from which they were submitted. Use root_dir instead.
        lines.append("set -x")
        lines.append("date '+%s'")
        lines.append("cd %s" % self.root_dir)
        # lines.append("module load fastp_0.20.1 samtools_1.12 minimap2_2.18")

        cmd = []
        cmd.append(self.fpmmp_path)
        cmd.append('-d %s' % self.root_dir)
        cmd.append('-D /Data/Fastq')
        cmd.append('-S {}\${SLURM_ARRAY_TASK_ID}'.format(self.trim_file))
        cmd.append('-p %s' % project_name)
        cmd.append('-C %s' % chemistry)
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
        lines.append("echo $?")
        lines.append("date '+%s'")

        with open(self.job_script_path, 'w') as f:
            logging.debug("Writing job script to %s" % self.job_script_path)
            for line in lines:
                # remove long spaces in some lines.
                f.write("%s\n" % line)


if __name__ == '__main__':
    import json

    hf2job = HumanFilter2Job('.', './good-sample-sheet.csv', '.', './fpmmp.sh', 16)
    results = hf2job._process_sample_sheet('./good-sample-sheet.csv')
    for key in results:
        print(key)
        print(json.dumps(results[key], indent=2))

'''
FIND=$(platform FIND)
  ### check filesystem for new sequence directory
  ### if exists, check for RTAComplete.txt files
pushd $seqpath/
path_count=(echo ${#raw_sequencing_loc[@]})
declare -a newdirs=$(echo "("; $FIND $seqpath/ -maxdepth 2 ! -path $seqpath/ -type d -mtime -1 ; echo ")")
count=${#newdirs[@]}
popd
labadmin_run=""
###sleep 5
echo sleeping
for dirpath in "${newdirs[@]}"
  do
    echo dirpath==$dirpath
		echo myfind==$myfind
		echo dirpath==$dirpath
    seqloc=$($myfind $dirpath -maxdepth 1 ! -path $dirpath -type f -name "RTAComplete.txt")

		#if [[ -e $dirpath/RTAComplete.txt ]]; then
    seqdir=$dirpath     ###$(dirname $seqloc)
		echo seqdir==${seqdir}
    fastq_output=$seqdir
    if [[ $seqloc ]]; then
      prep_data_loc "${seqdir}"
      #if [[ $filter_proc == "false" ]]; then
        if [[ $? == 0 ]]; then
          ###sleep 15
          process_data "${seqdir}" "$dirpath" ###"$curl"
        fi
      #fi
        #bcl_exit_status="0"
        if [[ $bcl_exit_status == "0" ]]; then
				  human_filter "${seqdir}" "$output_dir" "${fastq_output}"

          #fastqc_process "@"
        fi
    fi
		### end checking for RTAComplete.txt
  done

'''
