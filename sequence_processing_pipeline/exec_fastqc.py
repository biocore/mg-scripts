from multiprocessing import current_process
from subprocess import Popen, PIPE
import logging


def exec_fastqc(params):
    number_of_threads, file_path, output_path = params

    logging.debug(f"{current_process().name}: {number_of_threads} {file_path} "
                  f"{output_path}")

    # fastqc accepts n input files as parameters. However, legacy script
    # wrapped fastqc in a parallels process, and thus each process was
    # likely limited to a single input file.

    cmd = ['fastqc', '--noextract', '-t', number_of_threads, file_path, '-o',
           output_path]

    proc = Popen(' '.join(cmd), universal_newlines=True, shell=True,
                 stdout=PIPE, stderr=PIPE)

    stdout, stderr = proc.communicate()

    logging.debug(stdout)
    logging.debug(stderr)
    logging.debug(proc.returncode)

    # for now, return the file processed.
    return file_path
