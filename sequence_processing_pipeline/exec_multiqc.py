from subprocess import Popen, PIPE
import logging
from os import makedirs
from os.path import join
from functools import partial


def exec_multiqc(params):
    (number_of_threads, file_path, fastqc_output_dir, project_name,
     filter_type, fastq_trimmed_dir) = params

    multiqc_configuration_file_path = join('sequence_processing_pipeline',
                                           'multiqc-bclconvert-config.yaml')

    reports_results_dir = join(fastqc_output_dir, 'Reports')
    makedirs(reports_results_dir, exist_ok=True)

    base_path = partial(join, fastqc_output_dir, project_name)

    filtered_results_dir = base_path(filter_type)
    makedirs(filtered_results_dir, exist_ok=True)

    convert_results_dir = base_path('bclconvert')
    makedirs(convert_results_dir, exist_ok=True)

    multiqc_results_dir = base_path('multiqc')
    makedirs(multiqc_results_dir, exist_ok=True)

    raw_fastq_json_dir = base_path('json')

    json_dir = join(fastq_trimmed_dir, 'json')

    cmd = ['multiqc', '-c', multiqc_configuration_file_path, '--fullnames',
           '--force', reports_results_dir, raw_fastq_json_dir, json_dir,
           filtered_results_dir, convert_results_dir, '-o',
           multiqc_results_dir, '--interactive']

    cmd = ' '.join(cmd)

    proc = Popen(cmd, universal_newlines=True, shell=True, stdout=PIPE,
                 stderr=PIPE)

    stdout, stderr = proc.communicate()
    logging.debug(stdout)
    logging.debug(stderr)
    logging.debug(proc.returncode)
