import glob
import pgzip
import os
from sequence_processing_pipeline.util import iter_paired_files


def split_similar_size_bins(data_location_path, max_file_list_size_in_gb,
                            batch_prefix):
    '''Partitions input fastqs to coarse bins

    :param data_location_path: Path to the ConvertJob directory.
    :param max_file_list_size_in_gb: Upper threshold for file-size.
    :param batch_prefix: Path + file-name prefix for output-files.
    :return: The number of output-files created.
    '''

    # SPP root paths are of the form:
    # ../72dc6080-c453-418e-8a47-1ccd6d646a30/ConvertJob, and contain only
    # directories named after projects e.g:
    # 72dc6080-c453-418e-8a47-1ccd6d646a30/ConvertJob/Tell_Seq_15196.
    # Each of these directories contain R1 and R2 fastq files. Hence, path
    # is now the following:
    # add one more level to account for project_names nested under ConvertJob
    # dir.
    fastq_paths = glob.glob(data_location_path + '*/*/*.fastq.gz')

    # convert from GB and halve as we sum R1
    max_size = (int(max_file_list_size_in_gb) * (2 ** 30) / 2)

    split_offset = 0

    # ensure we are max-sized to start.
    current_size = max_size * 10
    fp = None

    for a, b in iter_paired_files(fastq_paths):
        r1_size = os.stat(a).st_size

        output_base = os.path.dirname(a).split('/')[-1]
        if current_size + r1_size > max_size:
            if fp is not None:
                fp.close()
            split_offset += 1
            current_size = r1_size
            fp = open(batch_prefix + '-%d' % split_offset, 'w')
        else:
            current_size += r1_size

        fp.write("%s\t%s\t%s\n" % (a, b, output_base))

    if fp is not None:
        fp.close()

    if split_offset == 0:
        raise ValueError("No splits made")

    return split_offset


def demux(id_map, fp, out_d, encoded, threads):
    """Split infile data based in provided map"""
    delimiter = '::MUX::'
    mode = 'wt'
    sep = '/'
    rec = '@'

    # load mapping information
    fpmap = {}
    for tmp in id_map:
        idx, r1, r2, outbase = tmp.strip().split('\t')
        fpmap[rec + idx] = (r1, r2, outbase)

    # this is our encoded file, and the one we care about
    idx = rec + encoded
    fname_r1, fname_r2, outbase = fpmap[idx]

    # setup output locations
    outdir = out_d + sep + outbase
    fullname_r1 = outdir + sep + fname_r1
    fullname_r2 = outdir + sep + fname_r2

    os.makedirs(outdir, exist_ok=True)
    current_fp_r1 = pgzip.open(fullname_r1, mode, thread=threads,
                               blocksize=2 * 10 ** 8)
    current_fp_r2 = pgzip.open(fullname_r2, mode, thread=threads,
                               blocksize=2 * 10 ** 8)

    current_fp = (current_fp_r1, current_fp_r2)

    # we assume R1 comes first... which is safe if data are coming
    # from samtools
    orientation_offset = 0

    # setup a parser
    id_ = iter(fp)
    seq = iter(fp)
    dumb = iter(fp)
    qual = iter(fp)

    for i, s, d, q in zip(id_, seq, dumb, qual):
        fname_encoded, id_ = i.split(delimiter, 1)

        if fname_encoded == idx:
            id_ = rec + id_

            current_fp[orientation_offset].write(id_)
            current_fp[orientation_offset].write(s)
            current_fp[orientation_offset].write(d)
            current_fp[orientation_offset].write(q)
            orientation_offset = 1 - orientation_offset  # switch between 0 / 1

    for outfp in current_fp:
        outfp.close()
