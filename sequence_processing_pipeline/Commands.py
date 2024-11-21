import glob
import gzip
import os
from sequence_processing_pipeline.util import iter_paired_files


def split_similar_size_bins(data_location_path, max_file_list_size_in_gb,
                            batch_prefix):
    '''Partitions input fastqs to coarse bins

    :param data_location_path: Path to the ConvertJob directory.
    :param max_file_list_size_in_gb: Upper threshold for file-size.
    :param batch_prefix: Path + file-name prefix for output-files.
    :return: The number of output-files created, size of largest bin.
    '''

    # SPP root paths are of the form:
    # ../72dc6080-c453-418e-8a47-1ccd6d646a30/ConvertJob, and contain only
    # directories named after projects e.g:
    # 72dc6080-c453-418e-8a47-1ccd6d646a30/ConvertJob/Tell_Seq_15196.
    # Each of these directories contain R1 and R2 fastq files. Hence, path
    # is now the following:
    # add one more level to account for project_names nested under ConvertJob
    # dir.
    # this will ignore the _I1_ reads that appear in the integrated result.
    fastq_paths = glob.glob(data_location_path + '/*/*_R?_001.fastq.gz')

    # convert from GB and halve as we sum R1
    max_size = (int(max_file_list_size_in_gb) * (2 ** 30) / 2)

    split_offset = 0

    # ensure we are max-sized to start.
    current_size = max_size * 10
    fp = None

    bucket_size = 0
    max_bucket_size = 0

    for a, b in iter_paired_files(fastq_paths):
        r1_size = os.stat(a).st_size
        r2_size = os.stat(b).st_size

        output_base = os.path.dirname(a).split('/')[-1]
        if current_size + r1_size > max_size:
            # bucket is full.
            if bucket_size > max_bucket_size:
                max_bucket_size = bucket_size

            # reset bucket_size.
            bucket_size = r1_size + r2_size

            if fp is not None:
                fp.close()

            split_offset += 1
            current_size = r1_size
            fp = open(batch_prefix + '-%d' % split_offset, 'w')
        else:
            # add to bucket_size
            bucket_size += r1_size + r2_size
            current_size += r1_size

        fp.write("%s\t%s\t%s\n" % (a, b, output_base))

    if fp is not None:
        fp.close()

    if split_offset == 0:
        raise ValueError("No splits made")

    return split_offset, max_bucket_size


def demux_cmd(id_map_fp, fp_fp, out_d, task, maxtask):
    with open(id_map_fp, 'r') as f:
        id_map = f.readlines()
        id_map = [line.strip().split('\t') for line in id_map]

    # fp needs to be an open file handle.
    # ensure task and maxtask are proper ints when coming from cmd-line.
    with open(fp_fp, 'r') as fp:
        demux(id_map, fp, out_d, int(task), int(maxtask))


def demux(id_map, fp, out_d, task, maxtask):
    """Split infile data based in provided map"""
    delimiter = '::MUX::'
    mode = 'wt'
    ext = '_001.fastq.gz'
    sep = '/'
    rec = '@'

    openfps = {}

    for offset, (idx, r1, r2, outbase) in enumerate(id_map):
        if offset % maxtask == task:
            idx = rec + idx

            # setup output locations
            outdir = out_d + sep + outbase
            fullname_r1 = outdir + sep + r1 + ext
            fullname_r2 = outdir + sep + r2 + ext

            os.makedirs(outdir, exist_ok=True)
            current_fp_r1 = gzip.open(fullname_r1, mode)
            current_fp_r2 = gzip.open(fullname_r2, mode)
            current_fp = {'1': current_fp_r1, '2': current_fp_r2}
            openfps[idx] = current_fp

    # setup a parser
    id_ = iter(fp)
    seq = iter(fp)
    dumb = iter(fp)
    qual = iter(fp)

    for i, s, d, q in zip(id_, seq, dumb, qual):
        fname_encoded, id_ = i.split(delimiter, 1)

        if fname_encoded not in openfps:
            continue

        orientation = id_[-2]  # -1 is \n
        current_fp = openfps[fname_encoded]
        id_ = rec + id_
        current_fp[orientation].write(id_)
        current_fp[orientation].write(s)
        current_fp[orientation].write(d)
        current_fp[orientation].write(q)

    for d in openfps.values():
        for f in d.values():
            f.close()
