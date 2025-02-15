#!/usr/bin/env python

import click
import os
import re
from os.path import abspath
from metapool.mp_strings import get_short_name_and_id
import re
import os
import gzip
import warnings
import pandas as pd
from glob import glob
from datetime import datetime
from metapool.mp_strings import get_short_name_and_id
from metapool.plate import PlateReplication
from collections import Counter
from string import ascii_letters, digits
from os import sep
from os.path import join, split, abspath
from collections import defaultdict

from metapool import (load_sample_sheet,
                      sample_sheet_to_dataframe, run_counts,
                      get_qiita_id_from_project_name)

def parse_illumina_run_id(run_id):
    """Parse a run identifier

    Parameters
    ----------
    run_id: str
        The name of a run

    Returns
    -------
    str:
        When the run happened (YYYY-MM-DD)
    str:
        Instrument code
    """

    # Format should be YYMMDD_machinename_XXXX_FC
    # this regex has two groups, the first one is the date, and the second one
    # is the machine name + suffix. This URL shows some examples
    # tinyurl.com/rmy67kw
    matches6 = re.search(r'^(\d{6})_(\w*)', run_id)

    # iSeq uses YYYYMMDD and has a trailing -XXXX value, for example
    # 20220303_FS10001773_6_BRB11606-1914
    # So the regex is updated to allow 8 numbers for the date, and to handle
    # the trailing -XXXX piece
    matches8 = re.search(r'^(\d{8})_([\w-]*)', run_id)

    if matches6 is None and matches8 is None:
        raise ValueError('Unrecognized run identifier format "%s". The '
                         'expected format is either '
                         'YYMMDD_machinename_XXXX_FC or '
                         'YYYYMMDD_machinename_XXXX-XXXX' %
                         run_id)

    if matches6 is None:
        matches = matches8
        fmt = "%Y%m%d"
    else:
        matches = matches6
        fmt = "%y%m%d"

    if len(matches.groups()) != 2:
        raise ValueError('Unrecognized run identifier format "%s". The '
                         'expected format is either '
                         'YYMMDD_machinename_XXXX_FC or '
                         'YYYYMMDD_machinename_XXXX-XXXX' %
                         run_id)

    # convert illumina's format to qiita's format
    run_date = datetime.strptime(matches[1], fmt).strftime('%Y-%m-%d')

    return run_date, matches[2]


# put together by Gail, based on the instruments we know of
INSTRUMENT_LOOKUP = pd.DataFrame({
    'FS10001773': {'machine prefix': 'FS', 'Vocab': 'Illumina iSeq',
                   'Machine type': 'iSeq', 'run_center': 'KLM'},
    'A00953': {'machine prefix': 'A', 'Vocab': 'Illumina NovaSeq 6000',
               'Machine type': 'NovaSeq', 'run_center': 'IGM'},
    'A00169': {'machine prefix': 'A', 'Vocab': 'Illumina NovaSeq 6000',
               'Machine type': 'NovaSeq', 'run_center': 'LJI'},
    'M05314': {'machine prefix': 'M', 'Vocab': 'Illumina MiSeq',
               'Machine type': 'MiSeq', 'run_center': 'KLM'},
    'K00180': {'machine prefix': 'K', 'Vocab': 'Illumina HiSeq 4000',
               'Machine type': 'HiSeq', 'run_center': 'IGM'},
    'D00611': {'machine prefix': 'D', 'Vocab': 'Illumina HiSeq 2500',
               'Machine type': 'HiSeq/RR', 'run_center': 'IGM'},
    'LH00444': {'machine prefix': 'LH', 'Vocab': 'Illumina NovaSeq X',
                'Machine type': 'NovaSeq', 'run_center': 'IGM'},
    'MN01225': {'machine prefix': 'MN', 'Vocab': 'Illumina MiniSeq',
                'Machine type': 'MiniSeq', 'run_center': 'CMI'}}).T


def get_machine_code(instrument_model):
    """Get the machine code for an instrument's code

    Parameters
    ----------
    instrument_model: str
        An instrument's model of the form A999999 or AA999999

    Returns
    -------
    """
    # the machine code represents the first 1 to 2 letters of the
    # instrument model
    machine_code = re.compile(r'^([a-zA-Z]{1,2})')
    matches = re.search(machine_code, instrument_model)
    if matches is None:
        raise ValueError('Cannot find a machine code. This instrument '
                         'model is malformed %s. The machine code is a '
                         'one or two character prefix.' % instrument_model)
    return matches[0]





def get_model_and_center(instrument_code):
    """Determine instrument model and center based on a lookup

    Parameters
    ----------
    instrument_code: str
        Instrument code from a run identifier.

    Returns
    -------
    str
        Instrument model.
    str
        Run center based on the machine's id.
    """
    run_center = "UCSDMI"
    instrument_model = instrument_code.split('_')[0]

    if instrument_model in INSTRUMENT_LOOKUP.index:
        run_center = INSTRUMENT_LOOKUP.loc[instrument_model, 'run_center']
        instrument_model = INSTRUMENT_LOOKUP.loc[instrument_model, 'Vocab']
    else:
        instrument_prefix = get_machine_code(instrument_model)

        if instrument_prefix not in INSTRUMENT_LOOKUP['machine prefix']:
            ValueError('Unrecognized machine prefix %s' % instrument_prefix)

        instrument_model = INSTRUMENT_LOOKUP[
            INSTRUMENT_LOOKUP['machine prefix'] == instrument_prefix
            ]['Vocab'].unique()[0]

    return instrument_model, run_center



def _find_filtered_files(fp):
    # assume that only fastq.gz files are desired and it's okay to ignore
    # other odd files that we find.

    # assume that all desired fastq.gz files will be stored in a directory
    # named 'filtered_sequences'. This will automatically ignore files in
    # 'only-adapter-filtered', 'zerofiles' and possibly other directories
    # in the future as well.

    # we expect the filepaths to be of the form:
    # <project_name_w_qiita_id>/filtered_sequences/<fastq_file>
    # However, we intentionally search across all directories to locate
    # any fastq files w/in a filtered_sequences subdirectory, even if
    # they don't match our expectation so that no potentially desirable
    # file is silently filtered out.
    files = glob(f"{fp}/**/filtered_sequences/*.fastq.gz", recursive=True)

    by_project = defaultdict(list)

    for fastq_fp in files:
        # generate the full path to each file for reference.
        full_path = abspath(fastq_fp)

        # process the relative path of each file to ensure each is of the
        # the expected form. Report anything out of the ordinary so that no
        # files are silently missed.

        # remove fp from the beginning of the file's path. The remaining path
        # should be of the form:
        # <project_name_w_qiita_id>/filtered_sequences/<fastq_file>
        tmp = fastq_fp.replace(fp, '')
        # remove any leading and/or trailing '/' characters from the
        # remaining path.
        # use os.sep instead of '/' to be more platform independent.
        tmp = tmp.strip(sep)
        tmp = tmp.split(sep)

        if len(tmp) != 3 or tmp[1] != 'filtered_sequences':
            raise ValueError(f"{fastq_fp} appears to be stored in an "
                             "unexpected location")

        # where tmp[0] is the project_name, keep a list of all associated with
        # that project.
        by_project[tmp[0]].append(full_path)

    return by_project




def _map_files_to_sample_ids(sample_ids, sample_files):
    # create a mapping of sample_ids to fastq filepaths.
    # organize processing by file to ensure each file is mapped and errors are
    # reported.
    results = defaultdict(list)

    for sample_file in sample_files:
        # make no assumptions about the naming format of fastq files, other
        # than their names will begin with a valid sample_id. This is to
        # support non-Illumina generated fastq files.
        _, file_name = split(sample_file)

        # for each sample_file, find the subset of sample_ids that match the
        # beginning of the file's name and assume that the longest matching
        # sample_id is the best and proper match. Ex: myfile_sre.fastq.gz
        # should match both the sample_ids [myfile_sre, myfile] but we want to
        # ensure myfile_sre is the correct match, and myfile_sre.fastq.gz is
        # not associated with myfile and myfile.fastq.gz is not associated with
        # myfile_sre.fastq.gz.
        matches = [s_id for s_id in sample_ids if file_name.startswith(s_id)]

        # sort the matches so the longest matching sample_id will be first.
        matches = sorted(matches, key=len, reverse=True)

        if not matches:
            # if no sample_ids could be matched to this file, then something
            # is wrong.
            raise ValueError(f"{sample_file} could not be matched to any "
                             "sample-id")

        # map sample_ids to a list of matching sample_files. We should expect
        # two matches for R1 and R2 files.
        results[matches[0]].append(sample_file)

    unmatched_sample_ids = set(sample_ids) - set(list(results.keys()))

    # after all possible sample_files have been matched to sample_ids, ensure
    # that the number of matching files for each sample_id is 2 and only 2,
    # assuming that one is an R1 file and the other is an R2 file.
    msgs = []
    for sample_id in results:
        count = len(results[sample_id])
        if count != 2:
            msgs.append(f"{sample_id} matched an unexpected number of samples "
                        f"({results[sample_id]})")
    if msgs:
        raise ValueError("\n".join(msgs))

    return results, unmatched_sample_ids


def _get_run_prefix(file_name):
    # aka forward, reverse, and indexed reads
    orientations = ['R1', 'R2', 'I1', 'I2']

    results = []

    # assume orientation is always present in the file's name.
    # assume that it is of one of the four forms above.
    # assume that it is always the right-most occurance of the four
    # orientations above.
    # assume that orientation is encapsulated with either '_' or '.'
    # e.g.: '_R1_', '.I2.'.
    # assume users can and will include any or all of the four
    # orientation as part of their filenames as well. e.g.:
    # ABC_7_04_1776_R1_SRE_S3_L007_R2_001.trimmed.fastq.gz
    for o in orientations:
        variations = [f"_{o}_", f".{o}."]
        for v in variations:
            # rfind searches from the end of the string, rather than
            # its beginning. It returns the position in the string
            # where the substring begins.
            results.append((file_name.rfind(v), o))

    # the orientation will be the substring found with the maximum
    # found value for pos. That is, it will be the substring that
    # begins at the rightest most position in the file name.
    results.sort(reverse=True)

    pos, orientation = results[0]

    # if no orientations were found, then return None.
    return None if pos == -1 else file_name[0:pos]




def process_sample(sample, prep_columns, run_center, run_date, run_prefix,
                   project_name, instrument_model, run_id, lane):
    # initialize result
    result = {c: '' for c in prep_columns}

    # hard-coded columns
    result["platform"] = "Illumina"
    result["sequencing_meth"] = "sequencing by synthesis"
    result["center_name"] = "UCSD"

    # manually generated columns
    result["run_center"] = run_center
    result["run_date"] = run_date
    result["run_prefix"] = run_prefix
    result["center_project_name"] = project_name
    result["instrument_model"] = instrument_model
    result["runid"] = run_id

    # lane is extracted from sample-sheet but unlike the others is passed
    # to this function explicitly.
    result["lane"] = lane

    # handle multiple types of sample-sheets, where columns such
    # as 'syndna_pool_number' may or may not be present.
    additional_columns = ['syndna_pool_number', 'mass_syndna_input_ng',
                          'extracted_gdna_concentration_ng_ul',
                          'vol_extracted_elution_ul',
                          'total_rna_concentration_ng_ul',
                          "sample_name", "experiment_design_description",
                          "library_construction_protocol", "sample_plate",
                          "i7_index_id", "index", "i5_index_id", "index2",
                          "sample_project", 'well_id_384', 'sample_well'
                          ]

    for attribute in additional_columns:
        if attribute in sample:
            result[attribute] = sample[attribute]

    # sanity-checks
    if 'well_id_384' in sample and 'sample_well' in sample:
        # sample_well will be 'Sample_Well' to the user.
        raise ValueError("'well_id_384' and 'Sample_Well' are both defined"
                         "in sample.")

    if 'well_id_384' not in sample and 'sample_well' not in sample:
        raise ValueError("'well_id_384' and 'Sample_Well' are both undefined"
                         "in sample.")

    well_id = result['well_id_384'] if 'well_id_384' in sample else result[
        'sample_well']

    result["well_description"] = '%s.%s.%s' % (sample.sample_plate,
                                               sample.sample_name,
                                               well_id)

    # return unordered result. let caller re-order columns as they see fit.
    return result



def agp_transform(frame, study_id):
    """If the prep belongs to the American Gut Project fill in some blanks

    Parameters
    ----------
    frame: pd.DataFrame
        The preparation file for a single project.
    study_id: str
        The Qiita study identifier for this preparations.

    Returns
    -------
    pd.DataFrame:
        If the study_id is "10317" then:
            - `center_name`, 'library_construction_protocol', and
              `experiment_design_description` columns are filled in with
              default values.
            - `sample_name` will be zero-filled no 9 digits.
        Otherwise no changes are made to `frame`.
    """
    if study_id == '10317':
        def zero_fill(name):
            if 'blank' not in name.lower() and name[0].isdigit():
                return name.zfill(9)
            return name

        frame['sample_name'] = frame['sample_name'].apply(zero_fill)
        frame['center_name'] = 'UCSDMI'
        frame['library_construction_protocol'] = 'Knight Lab KHP'
        frame['experiment_design_description'] = (
            'samples of skin, saliva and feces and other samples from the AGP')
    return frame


def preparations_for_run(run_id, path_to_fastq_files, sheet,
                          generated_prep_columns, carried_prep_columns):
    """Given a run's path and sample sheet generates preparation files

    Parameters
    ----------
    run_id: str
        Often times the name of the directory containing raw Illumina output.
    path_to_fastq_files: str
        Path to a directory of filtered fastq files, organized by project.
    sheet: dataFrame
        dataFrame() of SampleSheet contents
    generated_prep_columns: list
        List of required columns for output that are not expected in
        KLSampleSheet. It is expected that prep.py can derive these values
        appropriately.
    carried_prep_columns: list
        List of required columns for output that are expected in KLSampleSheet.
        Varies w/different versions of KLSampleSheet.

    Returns
    -------
    dict
        Dictionary keyed by run identifier, project name and lane. Values are
        preparations represented as DataFrames.
    """
    run_date, instrument_code = parse_illumina_run_id(run_id)
    instrument_model, run_center = get_model_and_center(instrument_code)

    fastq_by_project = _find_filtered_files(path_to_fastq_files)

    output = {}

    # well_description is no longer a required column, since the sample-name
    #  is now taken from column sample_name, which is distinct from sample_id.
    # sample_id remains the version of sample_name acceptable to bcl2fastq/
    #  bcl-convert, while sample_name is now the original string. The
    #  well_description and description columns are no longer needed or
    #  required, but a warning will be generated if neither are present as
    #  they are customarily present.

    # if 'well_description' is defined instead as 'description', rename it.
    # well_description is a recommended column but is not required.
    if 'well_description' not in set(sheet.columns):
        warnings.warn("'well_description' is not present in sample-sheet. It "
                      "is not a required column but it is a recommended one.")
        if 'description' in sheet:
            warnings.warn("Using 'description' instead of 'well_description'"
                          " because that column isn't present", UserWarning)
            # copy and drop the original column
            sheet['well_description'] = sheet['description'].copy()
            sheet.drop('description', axis=1, inplace=True)

    not_present = set(carried_prep_columns) - set(sheet.columns)

    if not_present:
        raise ValueError("Required columns are missing: %s" %
                         ', '.join(not_present))

    all_columns = sorted(carried_prep_columns + generated_prep_columns)

    for project, project_sheet in sheet.groupby('sample_project'):
        project_name, qiita_id = get_short_name_and_id(project)

        sample_files = fastq_by_project[project]
        sample_ids = sheet['sample_id'].tolist()

        results, unmapped_samples = _map_files_to_sample_ids(sample_ids,
                                                             sample_files)

        # make sure to put unmapped_samples into failed_samples.log

        # if the Qiita ID is not found then make for an easy find/replace
        if qiita_id == project:
            qiita_id = 'QIITA-ID'

        for lane, lane_sheet in project_sheet.groupby('lane'):
            # this is the portion of the loop that creates the prep
            data = []

            for well_id_col, sample in lane_sheet.iterrows():
                if isinstance(sample, pd.core.series.Series):
                    sample_id = well_id_col
                else:
                    sample_id = sample.sample_id

                if sample_id in unmapped_samples:
                    # ignore these as no matching fastq files were found
                    # for them, hence they didn't make it through conversion
                    # and host filtering.
                    continue

                if not sample_id in results:
                    raise ValueError("fastq.gz files were not found for "
                                     f"sample-id '{sample_id}'")

                fastq_files = results[sample_id]

                t1 = _get_run_prefix(fastq_files[0])
                t2 = _get_run_prefix(fastq_files[1])

                if t1 != t2:
                    raise ValueError(f"run_prefix values for '{sample_id}' "
                                     f"did not match: '{t1}', '{t2}'")

                run_prefix = t1

                data.append(process_sample(sample, all_columns, run_center,
                                           run_date, run_prefix, project_name,
                                           instrument_model, run_id, lane))

            if not data:
                warnings.warn('Project %s and Lane %s have no data' %
                              (project, lane), UserWarning)

            # the American Gut Project is a special case. We'll likely continue
            # to grow this study with more and more runs. So we fill some of
            # the blanks if we can verify the study id corresponds to the AGP.
            # This was a request by Daniel McDonald and Gail
            prep = agp_transform(pd.DataFrame(columns=all_columns,
                                              data=data), qiita_id)

            _check_invalid_names(prep.sample_name)

            output[(run_id, project, lane)] = prep

    return output

def _check_invalid_names(sample_names):
    # taken from qiita.qiita_db.metadata.util.get_invalid_sample_names
    valid = set(ascii_letters + digits + '.')

    def _has_invalid_chars(name):
        return bool(set(name) - valid)

    invalid = sample_names[sample_names.apply(_has_invalid_chars)]

    if len(invalid):
        warnings.warn('The following sample names have invalid '
                      'characters: %s' %
                      ', '.join(['"%s"' % i for i in invalid.values]))


@click.command()
@click.argument('run_dir', type=click.Path(exists=True, dir_okay=True,
                                           file_okay=False))
@click.argument('sample_sheet', type=click.Path(exists=True, dir_okay=False,
                                                file_okay=True))
@click.argument('output_dir', type=click.Path(writable=True))
@click.option('--verbose', help='list prep-file output paths, study_ids',
              is_flag=True)
def format_preparation_files(run_dir, sample_sheet, output_dir, verbose):
    """Generate the preparation files for the projects in a run

    RUN_DIR: should be the directory where the results of running bcl2fastq are
    saved.

    SAMPLE_SHEET: should be a CSV file that includes information for the
    samples and projects in RUN_DIR.

    OUTPUT_DIR: directory where the outputted preparations should be saved to.

    Preparations are stratified by project and by lane. Only samples with
    non-empty files are included. If "fastp-and-minimap2" is used, the script
    will collect sequence count stats for each sample and add them as columns
    in the preparation file.
    """
    sample_sheet = load_sample_sheet(sample_sheet)
    df_sheet = sample_sheet_to_dataframe(sample_sheet)

    stats = run_counts(run_dir, sample_sheet)
    stats['sample_name'] = \
        df_sheet.set_index('lane', append=True)['sample_name']

    # sample_sheet_to_dataframe() automatically lowercases the column names
    # before returning df_sheet. Hence, sample_sheet.CARRIED_PREP_COLUMNS also
    # needs to be lowercased for the purposes of tests in
    # preparation_for_run().
    c_prep_columns = [x.lower() for x in sample_sheet.CARRIED_PREP_COLUMNS]
    # returns a map of (run, project_name, lane) -> preparation frame
    preps = preparations_for_run(run_dir,
                                 df_sheet,
                                 sample_sheet.GENERATED_PREP_COLUMNS,
                                 c_prep_columns)

    os.makedirs(output_dir, exist_ok=True)

    for (run, project, lane), df in preps.items():
        fp = os.path.join(output_dir, f'{run}.{project}.{lane}.tsv')

        # stats are indexed by sample name and lane, lane is the first
        # level index. When merging, make sure to select the lane subset
        # that we care about, otherwise we'll end up with repeated rows
        df = df.merge(stats.xs(lane, level=1), how='left', on='sample_name')

        # strip qiita_id from project names in sample_project column
        df['sample_project'] = df['sample_project'].map(
            lambda x: re.sub(r'_\d+$', r'', x))

        # center_project_name is a legacy column that should mirror
        # the values for sample_project.
        df['center_project_name'] = df['sample_project']

        df.to_csv(fp, sep='\t', index=False)

        if verbose:
            # assume qiita_id is extractable and is an integer, given that
            # we have already passed error-checking.
            qiita_id = get_qiita_id_from_project_name(project)
            print("%s\t%s" % (qiita_id, abspath(fp)))

