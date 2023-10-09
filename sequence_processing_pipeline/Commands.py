import glob
import pgzip
import os
from os.path import split
import re


def split_similar_size_bins(data_location_path, max_file_list_size_in_gb,
                            batch_prefix):
    '''

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
    fastq_paths = glob.glob(data_location_path + '*/*.fastq.gz')
    fastq_paths.sort()

    # output_base_count is the number of fastq files in data_location_path.
    # fastq_paths should contain a list of R1 and R2 fastq files only.
    output_base_count = len(fastq_paths)

    files = iter(fastq_paths)

    # convert from GB and halve as we sum R1
    max_size = (int(max_file_list_size_in_gb) * (2 ** 30) / 2)

    r1 = iter(files)
    r2 = iter(files)

    split_offset = 0

    # ensure we are max-sized to start.
    current_size = max_size * 10
    fp = None

    r1_underscore = (re.compile(r'_R1_'), '_R1_', '_R2_')
    r1_dot = (re.compile(r'\.R1\.'), '.R1.', '.R2.')

    for a, b in zip(r1, r2):
        # a: /some_path/13722/115207/TCGA-05-4395-01A-01D-1203_110913_SN208_
        #  0238_BC008YACXX_s_8_rg.sorted.filtered.R1.trimmed.fastq.gz
        # b: /some_path/13722/115207/TCGA-05-4395-01A-01D-1203_110913_SN208_
        #  0238_BC008YACXX_s_8_rg.sorted.filtered.R2.trimmed.fastq.gz
        matched = False                                                             
        for pattern, r1_exp, r2_exp in (r1_underscore, r1_dot):                     
            if pattern.search(a):                                                  
                assert r2_exp in b                                                 
                assert a[:a.find(r1_exp)] == b[:b.find(r2_exp)]                 
                matched = True                                                      
                break                                                               
                                                                                
        if not matched:                                                             
            raise ValueError(f"Unable to match:\n{a}\n{b}") 

        r1_size = os.stat(a).st_size

        # output_base is the path up to the beginning of the filename.
        # we could do this in theory with os.path.split()[0] e.g.:
        # /panfs/dtmcdonald/human-depletion/t2t-only/13722/115207
        # /panfs/dtmcdonald/human-depletion/t2t-only/13722/115212
        output_base = split(a)

        if current_size + r1_size > max_size:
            if fp is not None:
                fp.close()
            split_offset += 1
            current_size = r1_size
            fp = open(batch_prefix + '-%d' % split_offset, 'w')
        else:
            current_size += r1_size

        fp.write("%s\t%s\t%s\n" % (a, b, output_base))
    fp.close()

    return split_offset


def demux_inline(argv1, argv2, argv3, argv4):
    delimiter = b'::MUX::'
    mode = 'wb'
    ext = b'.fastq.gz'
    sep = b'/'
    rec = b'@'

    # load mapping information
    fpmap = {}
    for l in open(argv1, 'rb'):
        idx, r1, r2, outbase = l.strip().split(b'\t')
        fpmap[rec + idx] = (
        r1, r2, outbase)  # append rec in mapping rather than per seq filter

    # gather other parameters
    fp = open(argv2, 'rb')
    out_d = argv3.encode('ascii')

    # this is our encoded file, and the one we care about
    idx = rec + argv4.encode('ascii')
    fname_r1, fname_r2, outbase = fpmap[idx]

    # setup output locations
    outdir = out_d + sep + outbase
    fullname_r1 = outdir + sep + fname_r1 + ext
    fullname_r2 = outdir + sep + fname_r2 + ext

    os.makedirs(outdir, exist_ok=True)
    current_fp_r1 = pgzip.open(fullname_r1, mode, thread=8,
                               blocksize=2 * 10 ** 8)
    current_fp_r2 = pgzip.open(fullname_r2, mode, thread=8,
                               blocksize=2 * 10 ** 8)

    current_fp = (current_fp_r1, current_fp_r2)

    # we assume R1 comes first...
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
