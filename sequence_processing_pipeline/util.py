import re


PAIR_UNDERSCORE = (re.compile(r'_R1_'), '_R1_', '_R2_')
PAIR_DOT = (re.compile(r'\.R1\.'), '.R1.', '.R2.')
PAIR_TESTS = (PAIR_UNDERSCORE, PAIR_DOT)


def iter_paired_files(files):
    """Yield matched r1/r2 paired files"""
    files = sorted(files)

    if len(files) % 2 != 0:
        raise ValueError("Files are not paired")

    files = iter(files)
    r1 = iter(files)
    r2 = iter(files)

    for r1_fp, r2_fp in zip(r1, r2):
        matched = False
        for pattern, r1_exp, r2_exp in PAIR_TESTS:
            if pattern.search(r1_fp):
                if r2_exp not in r2_fp:
                    raise ValueError(f"Cannot find '{r2_exp}' in '{r2_fp}'")

                # replace find w/find so that search for R1 and R2 begin
                # from the end of the string, not the beginning. This prevents
                # the code from breaking when filenames include R1 and R2 as
                # part of their name in addition to representing forward and
                # reversed reads e.g.:
                # LS_8_22_2014_R1_SRE_S3_L007_R1_001.trimmed.fastq.gz
                # LS_8_22_2014_R1_SRE_S3_L007_R2_001.trimmed.fastq.gz
                # using find(), r1_prefix and r2_prefix will be the following:
                # r1_prefix will be: LS_8_22_2014
                # r2_prefix will be: LS_8_22_2014_R1_SRE_S3_L007

                # r1_prefix = r1_fp[:r1_fp.find(r1_exp)]
                # r2_prefix = r2_fp[:r2_fp.find(r2_exp)]

                r1_prefix = r1_fp[:r1_fp.rfind(r1_exp)]
                r2_prefix = r2_fp[:r2_fp.rfind(r2_exp)]

                if r1_prefix != r2_prefix:
                    raise ValueError(f"Mismatch prefixes:\n{r1_prefix}\n"
                                     f"{r2_prefix}")

                matched = True
                break

        if not matched:
            raise ValueError(f"Unable to match:\n{r1_fp}\n{r2_fp}")

        yield (r1_fp, r2_fp)
