#!/usr/bin/env python

import click
import os
import pandas as pd
from os.path import abspath
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

from metapool import (
                      get_qiita_id_from_project_name)




PREP_MF_COLUMNS = ['sample_name', 'barcode', 'center_name',
                   'center_project_name', 'experiment_design_description',
                   'instrument_model', 'lane', 'library_construction_protocol',
                   'platform', 'run_center', 'run_date', 'run_prefix', 'runid',
                   'sample_plate', 'sequencing_meth', 'linker', 'primer',
                   'primer_plate', 'well_id_384', 'plating',
                   'extractionkit_lot', 'extraction_robot', 'tm1000_8_tool',
                   'primer_date', 'mastermix_lot', 'water_lot',
                   'processing_robot', 'tm300_8_tool', 'tm50_8_tool',
                   'project_name', 'orig_name', 'well_description',
                   'pcr_primers', 'target_gene', 'tm10_8_tool',
                   'target_subfragment', 'well_id_96']


REQUIRED_MF_COLUMNS = {'sample_name', 'barcode', 'primer', 'primer_plate',
                       'well_id_384', 'plating', 'extractionkit_lot',
                       'extraction_robot', 'tm1000_8_tool', 'primer_date',
                       'mastermix_lot', 'water_lot', 'processing_robot',
                       'tm300_8_tool', 'tm50_8_tool', 'sample_plate',
                       'project_name', 'orig_name', 'well_description',
                       'experiment_design_description', 'tm10_8_tool',
                       'library_construction_protocol', 'linker', 'platform',
                       'run_center', 'run_date', 'run_prefix', 'pcr_primers',
                       'sequencing_meth', 'target_gene', 'target_subfragment',
                       'center_name', 'center_project_name', 'well_id_96',
                       'instrument_model', 'runid'}

def get_run_prefix_mf(run_path, project):
    search_path = os.path.join(run_path, project, 'amplicon',
                               '*_SMPL1_S*R?_*.fastq.gz')
    results = glob(search_path)

    # at this stage there should only be two files forward and reverse
    if len(results) == 2:
        forward, reverse = sorted(results)
        if is_nonempty_gz_file(forward) and is_nonempty_gz_file(reverse):
            f, r = os.path.basename(forward), os.path.basename(reverse)
            if len(f) != len(r):
                raise ValueError("Forward and reverse sequences filenames "
                                 "don't match f:%s r:%s" % (f, r))

            # The first character that's different is the number in R1/R2. We
            # find this position this way because sometimes filenames are
            # written as _R1_.. or _R1.trimmed... and splitting on _R1 might
            # catch some substrings not part of R1/R2.
            for i in range(len(f)):
                if f[i] != r[i]:
                    i -= 2
                    break

            return f[:i]
        else:
            return None
    elif len(results) > 2:
        warnings.warn("%d possible matches found for determining run-prefix. "
                      "Only two matches should be found (forward and "
                      "reverse): %s" % (len(results),
                                        ', '.join(sorted(results))))

    return None



def preparations_for_run_mapping_file(run_path, mapping_file):
    """Given a run's path and mapping-file generates preparation files

        Parameters
        ----------
        run_path: str
            Path to the run folder
        mapping_file: pandas.DataFrame
            mapping-file to convert

        Returns
        -------
        dict
            Dictionary keyed by run identifier, project name and lane. Values
            are preparations represented as DataFrames.
        """

    # lowercase all columns
    mapping_file.columns = mapping_file.columns.str.lower()

    # add a faked value for 'lane' to preserve the original logic.
    # lane will always be '1' for amplicon runs.
    mapping_file['lane'] = pd.Series(
        ['1' for x in range(len(mapping_file.index))])

    # add a faked column for 'Sample_ID' to preserve the original logic.
    # count-related code will need to search run_directories based on
    # sample-id, not sample-name.
    mapping_file['Sample_ID'] = mapping_file.apply(
        # importing bcl_scrub_name led to a circular import
        lambda x: re.sub(r'[^0-9a-zA-Z\-\_]+', '_', x['sample_name']), axis=1)

    output = {}

    # lane is also technically a required column but since we generate it,
    # it's not included in REQUIRED_MF_COLUMNS.
    not_present = REQUIRED_MF_COLUMNS - set(mapping_file.columns)

    if not_present:
        raise ValueError("Required columns are missing: %s" %
                         ', '.join(not_present))

    for project, project_sheet in mapping_file.groupby('project_name'):
        _, qiita_id = get_short_name_and_id(project)

        # if the Qiita ID is not found then notify the user.
        if qiita_id == project:
            raise ValueError("Values in project_name must be appended with a"
                             " Qiita Study ID.")

        # note that run_prefix and run_id columns are required columns in
        # mapping-files. We expect these columns to be blank when seqpro is
        # run, however.
        run_prefix = get_run_prefix_mf(run_path, project)

        run_id = run_prefix.split('_SMPL1')[0]

        # return an Error if run_prefix could not be determined,
        # as it is vital for amplicon prep-info files. All projects will have
        # the same run_prefix.
        if run_prefix is None:
            raise ValueError("A run-prefix could not be determined.")

        for lane, lane_sheet in project_sheet.groupby('lane'):
            lane_sheet = lane_sheet.set_index('sample_name')

            # this is the portion of the loop that creates the prep
            data = []

            for sample_name, sample in lane_sheet.iterrows():
                row = {c: '' for c in PREP_MF_COLUMNS}

                row["sample_name"] = sample_name
                row["experiment_design_description"] = \
                    sample.experiment_design_description
                row["library_construction_protocol"] = \
                    sample.library_construction_protocol
                row["platform"] = sample.platform
                row["run_center"] = sample.run_center
                row["run_date"] = extract_run_date_from_run_id(run_id)
                row["run_prefix"] = run_prefix
                row["sequencing_meth"] = sample.sequencing_meth
                row["center_name"] = sample.center_name
                row["center_project_name"] = sample.center_project_name
                row["instrument_model"] = sample.instrument_model
                row["runid"] = run_id
                row["sample_plate"] = sample.sample_plate
                row["lane"] = lane
                row["barcode"] = sample.barcode
                row["linker"] = sample.linker
                row["primer"] = sample.primer
                row['primer_plate'] = sample.primer_plate
                if 'well_id_384' in sample:
                    row["well_id_384"] = sample.well_id_384
                elif 'sample_well' in sample:
                    row["sample_well"] = sample.sample_well
                row['well_id_96'] = sample.well_id_96
                row['plating'] = sample.plating
                row['extractionkit_lot'] = sample.extractionkit_lot
                row['extraction_robot'] = sample.extraction_robot
                row['tm1000_8_tool'] = sample.tm1000_8_tool
                row['primer_date'] = sample.primer_date
                row['mastermix_lot'] = sample.mastermix_lot
                row['water_lot'] = sample.water_lot
                row['processing_robot'] = sample.processing_robot
                row['tm300_8_tool'] = sample.tm300_8_tool
                row['tm50_8_tool'] = sample.tm50_8_tool
                row['tm10_8_tool'] = sample.tm10_8_tool
                row['project_name'] = sample.project_name
                row['orig_name'] = sample.orig_name
                row['well_description'] = sample.well_description
                row['pcr_primers'] = sample.pcr_primers
                row['target_gene'] = sample.target_gene
                row['target_subfragment'] = sample.target_subfragment
                data.append(row)

            if not data:
                warnings.warn('Project %s and Lane %s have no data' %
                              (project, lane), UserWarning)

            # the American Gut Project is a special case. We'll likely continue
            # to grow this study with more and more runs. So we fill some of
            # the blanks if we can verify the study id corresponds to the AGP.
            # This was a request by Daniel McDonald and Gail
            prep = agp_transform(pd.DataFrame(columns=PREP_MF_COLUMNS,
                                              data=data),
                                 qiita_id)

            _check_invalid_names(prep.sample_name)

            output[(run_id, project, lane)] = prep

    return output

@click.command()
@click.argument('run_dir', type=click.Path(exists=True, dir_okay=True,
                                           file_okay=False))
@click.argument('mapping_file', type=click.Path(exists=True, dir_okay=False,
                                                file_okay=True))
@click.argument('output_dir', type=click.Path(writable=True))
@click.option('--verbose', help='list prep-file output paths, study_ids',
              is_flag=True)
def format_preparation_files_mf(run_dir, mapping_file, output_dir, verbose):
    """Generate the preparation files for the projects in a run

    RUN_DIR: should be the directory where the results of running bcl2fastq are
    saved.

    MAPPING_FILE: should be a TSV file that includes information for the
    samples and projects in RUN_DIR.

    OUTPUT_DIR: directory where the outputted preparations should be saved to.

    Preparations are stratified by project and by lane. Only samples with
    non-empty files are included.
    """
    df_mf = pd.read_csv(mapping_file, delimiter='\t')

    # run_counts cannot be determined for a single multiplexed file.

    # returns a map of (run, project_name, lane) -> preparation frame
    preps = preparations_for_run_mapping_file(run_dir, df_mf)

    os.makedirs(output_dir, exist_ok=True)

    for (run, project, lane), df in preps.items():
        fp = os.path.join(output_dir, f'{run}.{project}.{lane}.tsv')

        df.to_csv(fp,
                  sep='\t',
                  index=False,
                  # finalize column order

                  columns=['sample_name', 'center_name', 'center_project_name',
                           'experiment_design_description', 'instrument_model',
                           'lane', 'library_construction_protocol',
                           'platform', 'run_center', 'run_date', 'run_prefix',
                           'runid', 'sample_plate', 'sequencing_meth',
                           'barcode', 'linker', 'primer', 'extraction_robot',
                           'extractionkit_lot', 'mastermix_lot', 'orig_name',
                           'pcr_primers', 'plating', 'primer_date',
                           'primer_plate', 'processing_robot', 'project_name',
                           'target_gene', 'target_subfragment',
                           'tm1000_8_tool', 'tm300_8_tool', 'tm50_8_tool',
                           'tm10_8_tool', 'water_lot', 'well_description',
                           'well_id_384', 'well_id_96'])

        if verbose:
            # assume qiita_id is extractable and is an integer, given that
            # we have already passed error-checking.
            qiita_id = get_qiita_id_from_project_name(project)
            print("%s\t%s" % (qiita_id, abspath(fp)))


if __name__ == '__main__':
    format_preparation_files_mf()

