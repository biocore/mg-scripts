from metapool import KLSampleSheet, validate_and_scrub_sample_sheet
from os import stat, listdir, makedirs
from os.path import join, basename
from sequence_processing_pipeline.Job import Job
from sequence_processing_pipeline.PipelineError import PipelineError
from sequence_processing_pipeline.Pipeline import Pipeline
from shutil import move
import logging
from json import dumps
from datetime import date
from sequence_processing_pipeline.Commands import split_similar_size_bins


logging.basicConfig(level=logging.DEBUG)


'''
TODO:
The only two methods not called in testing are:
in _filter
in run
'''


class NuQCJob(Job):
    def __init__(self, fastq_root_dir, output_path, sample_sheet_path,
                 minimap_database_paths, queue_name,
                 node_count, nprocs, wall_time_limit, jmem, fastp_path,
                 minimap2_path, samtools_path, modules_to_load, qiita_job_id,
                 pool_size, max_array_length):
        """
        Submit a Torque job where the contents of fastq_root_dir are processed
        using fastp, minimap, and samtools. Human-genome sequences will be
        filtered out, if needed.
        :param fastq_root_dir: Path to a dir of Fastq files, org. by project.
        :param output_path: Path where all pipeline-generated files live.
        :param sample_sheet_path: Path to a sample sheet file.
        :param minimap_database_paths: Path to human genome databases in env.
        :param queue_name: Torque queue name to use in running env.
        :param node_count: Number of nodes to use in running env.
        :param nprocs: Number of processes to use in runing env.
        :param wall_time_limit: Hard wall-clock-time limit (in min) for procs.
        :param jmem: String representing total memory limit for entire job.
        :param fastp_path: The path to the fastp executable
        :param minimap2_path: The path to the minimap2 executable
        :param samtools_path: The path to the samtools executable
        :param modules_to_load: A list of Linux module names to load
        :param qiita_job_id: identify Torque jobs using qiita_job_id
        :param pool_size: The number of jobs to process concurrently.

        """
        super().__init__(fastq_root_dir,
                         output_path,
                         'QCJob',
                         [fastp_path, minimap2_path, samtools_path],
                         max_array_length,
                         modules_to_load=modules_to_load)
        self.sample_sheet_path = sample_sheet_path
        self._file_check(self.sample_sheet_path)
        metadata = self._process_sample_sheet()
        self.sample_ids = metadata['sample_ids']
        self.project_data = metadata['projects']
        self.needs_trimming = metadata['needs_adapter_trimming']
        self.nprocs = 16 if nprocs > 16 else nprocs
        self.chemistry = metadata['chemistry']
        # instead of a list of individual .mmi files,  minimap_database_paths
        # is now a path to a directory of .mmi files.
        self.minimap_database_paths = minimap_database_paths
        self.queue_name = queue_name
        self.node_count = node_count
        self.wall_time_limit = wall_time_limit
        self.jmem = jmem
        self.fastp_path = fastp_path
        self.minimap2_path = minimap2_path
        self.samtools_path = samtools_path
        self.qiita_job_id = qiita_job_id
        self.pool_size = pool_size
        self.suffix = 'fastq.gz'

        ### BEGIN NEW PARAMS
        # sizes in GB
        self.max_file_list_size_in_gb = 2  # max size of GB data to process in a job
        self.data_location = "/panfs/dtmcdonald/human-depletion/t2t-only/*/*"

        datestamp = date.today().strftime("%Y.%m.%d")

        self.batch_prefix = join(self.output_path,
                                 f"1-hd-split-pangenome-pe-{datestamp}")

        self.output = "/panfs/dtmcdonald/human-depletion/pangenome-adapter-filter"
        self.tmpdir = "/scratch/tcga_compute/tmp"
        ### END NEW PARAMS

        self.minimum_bytes = 3100

        if not isinstance(self.needs_trimming, bool):
            raise ValueError("needs_adapter_trimming must be boolean.")

        # Validate project settings in [Bioinformatics]
        for project in self.project_data:
            if project['ForwardAdapter'] == 'NA':
                project['ForwardAdapter'] = None

            if project['ReverseAdapter'] == 'NA':
                project['ReverseAdapter'] = None

            if project['ForwardAdapter'] is None:
                if project['ReverseAdapter'] is not None:
                    raise ValueError(("ForwardAdapter is declared but not "
                                      "ReverseAdapter."))

            if project['ReverseAdapter'] is None:
                if project['ForwardAdapter'] is not None:
                    raise ValueError(("ReverseAdapter is declared but not "
                                      "ForwardAdapter."))

            if not isinstance(project['HumanFiltering'], bool):
                raise ValueError("needs_adapter_trimming must be boolean.")

    def _filter_empty_fastq_files(self, filtered_directory,
                                  empty_files_directory,
                                  minimum_bytes):
        '''
        Filters out and moves fastq files that are below threshold.
        :param filtered_directory:
        :param empty_files_directory:
        :param minimum_bytes:
        :return:
        '''
        empty_list = []

        for entry in listdir(filtered_directory):
            if '_R1_' in entry:
                reverse_entry = entry.replace('_R1_', '_R2_')
                full_path = join(filtered_directory, entry)
                full_path_reverse = join(filtered_directory, reverse_entry)
                if stat(full_path).st_size <= minimum_bytes or stat(
                        full_path_reverse).st_size <= minimum_bytes:
                    logging.debug(f'moving {entry} and {reverse_entry}'
                                  f' to empty list.')
                    empty_list.append(full_path)
                    empty_list.append(full_path_reverse)

        if empty_list:
            logging.debug(f'making directory {empty_files_directory}')
            makedirs(empty_files_directory, exist_ok=True)

        for item in empty_list:
            logging.debug(f'moving {item}')
            move(item, empty_files_directory)

    def run(self, callback=None):
        # now a single job-script will be created to process all projects at
        # the same time, and intelligently handle adapter-trimming as needed
        # as well as human-filtering.
        job_script_path = self._generate_job_script()

        batch_count = split_similar_size_bins(self.data_location,
                                              self.max_file_list_size_in_gb,
                                              self.batch_prefix)

        export_params = [f"MMI={self.minimap_database_paths}",
                         f"PREFIX={self.batch_prefix}",
                         f"OUTPUT={self.output}", f"TMPDIR={self.tmpdir}"]

        job_params = ['-J', self.batch_prefix, '--array', batch_count,
                      '--mem', '20G', '--export', ','.join(export_params)]

        # job_script_path formerly known as:
        #  process.multiprep.pangenome.adapter-filter.pe.sbatch
        job_info = self.submit_job(job_script_path,
                                   # job_parameters - was None
                                   ' '.join(job_params),
                                   # script_parameters
                                   None,
                                   # assume we want to exec from the log_path
                                   # for now.
                                   exec_from=self.log_path,
                                   callback=callback)

        job_id = job_info['job_id']
        logging.debug(f'QCJob {job_id} completed')

        '''
        for project in self.project_data:
            project_name = project['Sample_Project']
            needs_human_filtering = project['HumanFiltering']

            source_dir = join(self.output_path, project_name)


            # TODO: UNCONFIRMED COMPLETED FILES WILL BE CREATED
            # confirm Slurm job was successful by using .completed files.
            # if not self._get_failed_indexes(project_name, job_id):
            #     raise PipelineError("QCJob did not complete successfully.")

            ### DANIEL'S NEW CODE DOES NOT APPEAR TO FILTER FOR
                ZERO-LENGTH FILES ###

            # determine where the filtered fastq files can be found and move
            # the 'zero-length' files to another directory.
            if needs_human_filtering is True:
                filtered_directory = join(source_dir, 'filtered_sequences')
            else:
                filtered_directory = join(source_dir, 'trimmed_sequences')

            empty_files_directory = join(source_dir, 'zero_files')
            self._filter_empty_fastq_files(filtered_directory,
                                           empty_files_directory,
                                           self.minimum_bytes)
        '''

    def _process_sample_sheet(self):
        sheet = KLSampleSheet(self.sample_sheet_path)
        valid_sheet = validate_and_scrub_sample_sheet(sheet)

        if not valid_sheet:
            s = "Sample sheet %s is not valid." % self.sample_sheet_path
            raise PipelineError(s)

        header = valid_sheet.Header
        chemistry = header['chemistry']

        if header['Assay'] not in Pipeline.assay_types:
            s = "Assay value '%s' is not recognized." % header['Assay']
            raise PipelineError(s)

        needs_adapter_trimming = (
                    header['Assay'] in [Pipeline.METAGENOMIC_PTYPE,
                                        Pipeline.METATRANSCRIPTOMIC_PTYPE])

        sample_ids = []
        for sample in valid_sheet.samples:
            sample_ids.append((sample['Sample_ID'], sample['Sample_Project']))

        bioinformatics = valid_sheet.Bioinformatics

        # reorganize the data into a list of dictionaries, one for each row.
        # the ordering of the rows will be preserved in the order of the list.
        lst = bioinformatics.to_dict('records')

        # convert true/false and yes/no strings to true boolean values.
        for record in lst:
            for key in record:
                if record[key].strip().lower() in ['true', 'yes']:
                    record[key] = True
                elif record[key].strip().lower() in ['false', 'no']:
                    record[key] = False

        # human-filtering jobs are scoped by project. Each job requires
        # particular knowledge of the project.
        return {'chemistry': chemistry,
                'projects': lst,
                'sample_ids': sample_ids,
                'needs_adapter_trimming': needs_adapter_trimming
                }

    def _generate_job_script(self):
        # TODO: Move environment variables into params?
        lines = []
        lines.append("#!/bin/bash -l\n")
        lines.append("#SBATCH -J pangenome-filtering\n")
        lines.append("#SBATCH --time 96:00:00\n")
        lines.append("#SBATCH --mem 20gb\n")
        lines.append("#SBATCH -N 1\n")
        lines.append("#SBATCH -c 8\n")
        lines.append("#SBATCH --output %x-%A_%a.out\n")
        lines.append("#SBATCH --error %x-%A_%a.err\n")
        lines.append("\n")
        lines.append("if [[ -z \"${SLURM_ARRAY_TASK_ID}\" ]]; then\n")
        lines.append("\techo \"Not operating within an array\"\n")
        lines.append("\texit 1\n")
        lines.append("fi\n")
        lines.append("\n")
        lines.append("mamba activate human-depletion\n")
        lines.append("\n")
        lines.append("set -x\n")
        lines.append("set -e\n")
        lines.append("\n")
        lines.append("# set a temp directory, make a new unique one under it"
                     " and\n")
        lines.append("# make sure we clean up as we're dumping to shm\n")
        lines.append("# DO NOT do this casually. Only do a clean up like this"
                     " if\n")
        lines.append("# you know for sure TMPDIR is what you want.\n")
        lines.append("\n")
        lines.append("mkdir -p {TMPDIR}\n")
        lines.append("export TMPDIR=${TMPDIR}\n")
        lines.append("export TMPDIR=$(mktemp -d)\n")
        lines.append("echo $TMPDIR\n")
        lines.append("\n")
        lines.append("function cleanup {\n")
        lines.append("\techo \"Removing $TMPDIR\"\n")
        lines.append("\trm -fr $TMPDIR\n")
        lines.append("\tunset TMPDIR\n")
        lines.append("}\n")
        lines.append("trap cleanup EXIT\n")
        lines.append("\n")
        lines.append("export FILES=$(pwd)/$(printf \"%s-%d\" ${PREFIX} "
                     "${SLURM_ARRAY_TASK_ID})\n")
        lines.append("if [[ ! -f ${FILES} ]]; then\n")
        lines.append("\tlogger ${FILES} not found\n")
        lines.append("\texit 1\n")
        lines.append("fi\n")
        lines.append("\n")
        lines.append("delimiter=::MUX::\n")
        lines.append("n=$(wc -l ${FILES} | cut -f 1 -d\" \")\n")
        lines.append("\n")
        lines.append("for i in $(seq 1 ${n})\n")
        lines.append("do\n")
        lines.append("\tline=$(head -n ${i} ${FILES} | tail -n 1)\n")
        lines.append("\tr1=$(echo ${line} | cut -f 1 -d\" \")\n")
        lines.append("\tr2=$(echo ${line} | cut -f 2 -d\" \")\n")
        lines.append("\tbase=$(echo ${line} | cut -f 3 -d\" \")\n")
        lines.append("\tr1_name=$(basename ${r1} .fastq.gz)\n")
        lines.append("\tr2_name=$(basename ${r2} .fastq.gz)\n")
        lines.append("\n")
        lines.append("\techo \"${i}	${r1_name}	${r2_name}	${base}\" >> "
                     "${TMPDIR}/id_map\n")
        lines.append("\n")
        lines.append("\tfastp \\\n")
        lines.append("\t\t-l 45 \\\n")
        lines.append("\t\t-i ${r1} \\\n")
        lines.append("\t\t-I ${r2} \\\n")
        lines.append("\t\t-w 7 \\\n")
        lines.append("\t\t--adapter_fasta fastp_known_adapters_formatted.fna"
                     " \\\n")
        lines.append("\t\t--html /dev/null \\\n")
        lines.append("\t\t--json /dev/null \\\n")
        lines.append("\t\t--stdout | \\\n")
        lines.append("\t\tsed -r \"1~4s/^@(.*)/@${i}${delimiter}\\1/\"\n")
        lines.append("done > ${TMPDIR}/seqs.fastq\n")
        lines.append("\n")
        lines.append("for mmi in ${MMI}/*.mmi\n")
        lines.append("do\n")
        lines.append("\techo \"$(date) :: $(basename ${mmi})\"\n")
        lines.append("\tminimap2 -2 -ax sr -t 7 ${mmi} ${TMPDIR}/seqs.fastq |"
                     " \\\n")
        lines.append("\t\tsamtools fastq -@ 1 -f 12 -F 256 > "
                     "${TMPDIR}/seqs_new.fastq\n")
        lines.append("\techo $(du -sh ${TMPDIR})\n")
        lines.append("\tmv ${TMPDIR}/seqs_new.fastq ${TMPDIR}/seqs.fastq\n")
        lines.append("done\n")
        lines.append("\n")
        lines.append("mkdir -p ${OUTPUT}\n")
        lines.append("\n")
        lines.append("function runner () {\n")
        lines.append("\ti=${1}\n")
        lines.append("\tpython demux-inline.py ${TMPDIR}/id_map "
                     "${TMPDIR}/seqs.fastq ${OUTPUT} ${i}\n")
        lines.append("}\n")
        lines.append("\n")
        lines.append("export -f runner\n")
        lines.append("jobs=${SLURM_JOB_CPUS_PER_NODE}\n")
        lines.append("\n")
        lines.append("echo \"$(date) :: demux start\"\n")
        lines.append("# let it do its thing\n")
        lines.append("seq 1 ${n} | parallel -j ${jobs} runner\n")
        lines.append("echo \"$(date) :: demux stop\"\n")
        lines.append("\n")
        lines.append("rm -fr ${TMPDIR}\n")

        job_script_path = join(self.output_path, 'process_all_fastq_files.sh')

        with open(job_script_path, 'w') as f:
            f.write('\n'.join(lines))

        return job_script_path
