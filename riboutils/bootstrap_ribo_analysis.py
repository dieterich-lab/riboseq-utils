#! /usr/bin/env python3

import argparse
import os
import pandas as pd
import yaml

import misc.parallel as parallel
import misc.utils as utils

import logging
import misc.logging_utils as logging_utils
logger = logging.getLogger(__name__)

sample_sheet_file_type_choices = [
    'AUTO',
    'excel',
    'hdf5',
    'csv'
]
default_sample_sheet_file_type = 'AUTO'
default_sheet_name = "Sheet1"
default_riboseq_sample_types = ['riboseq']
default_rnaseq_sample_types = ['rna-seq', 'total-rna-seq']

def get_sample_name(condition:str, sample_type:str, cell_type:str=None, 
        replicate:str=None, lane:str=None):
    """ Construct a sample name by concatenating the respective properties of 
    the sample.
    
    The format of the name is:
        <condition>.<sample_type>[.cell-type-<cell_type>][.rep-<replicate>][.lane-<lane>]
        
    The optional parts are skipped if the respective value is None.
    
    Parameters
    ----------
    condition: string
        The name of the condition for the sample, e.g., "sham.cm"
        
    sample_type: string
        The type of the sample, e.g., "riboseq"
        
    cell_type: string
        The type of cell (tissue, etc.) from which the sample came, e.g., "cm"
        
    replicate: string
        An identifier for the (biological) replicate, e.g., "mouse-403"
        
    lane: string
        An identifier for the lane of the sample, e.g., "2"
        
    Returns
    -------
    sample_name: string, or None
        The name, constructed as indicated above. If condition is None, NaN, or
        a zero-length string, then None is returned.
    """
    if pd.isnull(condition):
        return None
    
    if condition is None:
        return None

    if len(condition) == 0:
        return None
    
    sample_name = "{}.{}".format(str(condition), str(sample_type))
    
    if cell_type is not None:
        sample_name = "{}.cell-type-{}".format(sample_name, str(cell_type))
    
    if replicate is not None:
        sample_name = "{}.rep-{}".format(sample_name, str(replicate))
        
    if lane is not None:
        sample_name = "{}.lane-{}".format(sample_name, str(lane))
        
    return sample_name

def _get_full_condition_name_helper(row):
    condition = row.get("condition")
    sample_type = row.get("sample_type")
    cell_type = row.get("cell_type")
    
    full_condition_name = get_sample_name(condition, sample_type, cell_type)
    return full_condition_name

def _get_replicate_name_helper(row):
    condition = row.get("condition")
    sample_type = row.get("sample_type")
    cell_type = row.get("cell_type")
    replicate = row.get("replicate")
    
    replicate_name = get_sample_name(
        condition, 
        sample_type, 
        cell_type, 
        replicate
    )

    return replicate_name

def _get_replicate_filename_helper(row):    
    location = row.get("location")
    replicate_name = row.get("replicate_name")
    
    if replicate_name is None:
        return None
    
    replicate_filename = "{}.fastq.gz".format(replicate_name)
    replicate_filename = os.path.join(location, replicate_filename)
    return replicate_filename

def _get_sample_name_helper(row):
    condition = row.get("condition")
    sample_type = row.get("sample_type")
    cell_type = row.get("cell_type")
    replicate = row.get("replicate")
    lane = row.get("lane")
    
    sample_name = get_sample_name(
        condition, 
        sample_type, 
        cell_type, 
        replicate, 
        lane
    )

    return sample_name

def _get_sample_filename_helper(row):    
    location = row.get("location")
    sample_name = row.get("sample_name")
    
    if sample_name is None:
        return None
    
    sample_filename = "{}.fastq.gz".format(sample_name)
    sample_filename = os.path.join(location, sample_filename)
    return sample_filename


def _create_symlink(row, overwrite):
    # first, ensure the relevant paths are valid
    if pd.isnull(row['filename']):
        return None

    if not os.path.exists(row['sample_filename']):
        msg = "Could not find original sample: {}".format(
            row['sample_filename'])
        raise FileNotFoundError(msg)

    # create the necessary directory structure, but do not overwrite any
    # existing file, since it could be "real" data
    utils.create_symlink(row['sample_filename'], row['filename'], 
        remove=overwrite, create=True)

def pool_lanes(g):
    to_cat = g['filename']
    row = g.iloc[0]
    rep_file = row['replicate_filename']
    utils.concatenate_files(to_cat, rep_file)
    
def get_condition_replicates(g):
    row = g.iloc[0]
    full_condition_name = row['full_condition_name']
    
    all_replicates = g['replicate_name'].unique()
    all_replicates = list(all_replicates)
    
    return {full_condition_name: all_replicates}

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Create symlinks and partial config file for samples and "
        "replicates from a csv sample sheet. The csv file must include the "
        "following columns: sample_filename, condition, sample_type. "
        "Optionally, it can also use the following columns in the filenames: "
        "cell_type, replicate_name, lane. The format of the symlinks and sample "
        "names is as follows:<condition>.<sample_type>[.cell-type-<cell_type>]"
        "[.rep-<replicate>][.lane-<lane>]. The optional parts are skipped if "
        "they are not present in the sample sheet. The script first "
        "concatenates samples with the same condition, sample type, cell type "
        "and replicate identifiers (but with different lanes). The "
        "\"biological_replicates\" group samples with the same condition, "
        "sample type and cell type.\n\nThe sample sheet can also contain "
        "additional columns, but they are ignored.")

    parser.add_argument('sample_sheet', help="The csv sample sheet. It can "
        "also be an excel file if it has the file extension \"xls\" or "
        "\"xlsx\" or hdf5 if the filetype is \"hdf\", \"hdf5\", \"h5\" or "
        "\"he5\".")

    parser.add_argument('out', help="The (partial) yaml config file created "
        "based on the sample sheet")

    parser.add_argument('--sample-sheet-file-type', help="The file type of "
        "the sample sheet. By default (\"AUTO\"), this is guessed based on "
        "the extension.", choices=sample_sheet_file_type_choices,
        default=default_sample_sheet_file_type)

    parser.add_argument('--sheet-name', help="For excel and hdf5 files, the "
        "name of the sheet in the workbook which contains the sample sheet.",
        default=default_sheet_name)

    parser.add_argument('--riboseq-sample-types', help="The \"sample_type\"s "
        "to treat as riboseq samples", nargs='*', 
        default=default_riboseq_sample_types)

    parser.add_argument('--rnaseq-sample-types', help="The \"sample_type\"s "
        "to treat as rna-seq samples", nargs='*',
        default=default_rnaseq_sample_types)

    parser.add_argument('--overwrite', help="If this flag is given, then "
        "files at the symlink locations will be overwritten. N.B. THIS COULD "
        "DESTROY THE ORIGINAL DATA FILES! BE CAREFUL!!!", action='store_true')

    parser.add_argument('--no-symlinks', help="If this flag is given, then "
        "symlinks will not be created. Namely, only the yaml config file will "
        "be written.", action='store_true')
    
    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)

    msg = "Reading sample sheet"
    logger.info(msg)

    # make sure we convert to the correct data type
    converters = {
        "condition": str,
        "sample_type": str,
        "cell_type": str,
        "replicate": str,
        "lane": str,
        "random_field": str
    }

    sample_sheet = utils.read_df(args.sample_sheet, 
        filetype=args.sample_sheet_file_type, 
        skip_blank_lines=True, 
        sheet=args.sheet_name, 
        converters=converters
    )

    # make sure we have the necessary columns
    if 'condition' not in sample_sheet.columns:
        msg = "\"condition\" is not present in the sample sheet"
        raise ValueError(msg)

    if 'sample_type' not in sample_sheet.columns:
        msg = "\"sample_type\" is not present in the sample sheet"
        raise ValueError(msg)

    if 'sample_filename' not in sample_sheet.columns:
        msg = "\"sample_filename\" is not present in the sample sheet"
        raise ValueError(msg)

    msg = "Creating filenames"
    logger.info(msg)

    sample_sheet['sample_name'] = parallel.apply_df_simple(
        sample_sheet, 
        _get_sample_name_helper
    )

    sample_sheet['filename'] = parallel.apply_df_simple(
        sample_sheet, 
        _get_sample_filename_helper
    )

    if not args.no_symlinks:
        msg = "Creating symlinks"
        logger.info(msg)

        parallel.apply_df_simple(sample_sheet, _create_symlink, args.overwrite)

    msg = "Pooling samples for each replicate from different lanes"
    logger.info(msg)

    sample_sheet['replicate_name'] = parallel.apply_df_simple(
        sample_sheet, 
        _get_replicate_name_helper
    )

    sample_sheet['replicate_filename'] = parallel.apply_df_simple(
        sample_sheet, 
        _get_replicate_filename_helper
    )
        
    replicate_groups = sample_sheet.groupby('replicate_name')

    if not args.no_symlinks:
        replicate_groups.apply(pool_lanes)

    msg = "Extracting replicate names for config"
    logger.info(msg)

    # finally, create the yaml config file
    m_riboseq = sample_sheet['sample_type'].isin(args.riboseq_sample_types)
    riboseq_samples = utils.dataframe_to_dict(
        sample_sheet[m_riboseq], 
        'replicate_name', 
        'replicate_filename'
    )

    m_rnaseq = sample_sheet['sample_type'].isin(args.rnaseq_sample_types)
    rnaseq_samples = utils.dataframe_to_dict(
        sample_sheet[m_rnaseq], 
        'replicate_name', 
        'replicate_filename'
    )

    msg = "Grouping replicates by condition"
    logger.info(msg)

    # so there is not a random magic number later...
    num_procs = 1

    sample_sheet['full_condition_name'] = parallel.apply_df_simple(
        sample_sheet, 
        _get_full_condition_name_helper
    )
        
    ribo_condition_groups = sample_sheet[m_riboseq].groupby('full_condition_name')
    ribo_condition_groups = parallel.apply_parallel_groups(
        ribo_condition_groups, 
        num_procs, 
        get_condition_replicates
    )

    ribo_condition_groups = utils.merge_dicts(*ribo_condition_groups)

    rna_condition_groups = sample_sheet[m_rnaseq].groupby('full_condition_name')
    rna_condition_groups = parallel.apply_parallel_groups(
        rna_condition_groups, 
        num_procs, 
        get_condition_replicates
    )

    rna_condition_groups = utils.merge_dicts(*rna_condition_groups)

    msg = "Writing partial config file"
    logger.info(msg)
    config = {
        'project_name': "<PROJECT_NAME>",
        'note': "<NOTE>",
        'gtf': "<GTF>",
        'fasta': "<FASTA>",
        'star_index': "<STAR_INDEX>",
        'ribosomal_index': "<RIBOSOMAL_INDEX>",
        'ribosomal_fasta': "<RIBOSOMAL_FASTA>",
        'genome_base_path': "<GENOME_BASE_PATH>",
        'genome_name': "<GENOME_NAME>",
        'orf_note': "<ORF_NOTE>",
        'adapter_file': "<RIBO_ADAPTER_FILE>",
        'rna_adapter_file': "<RNA_ADAPTER_FILE>",
        'riboseq_data': "<RIBO_DATA_PATH>",
        'rnaseq_data': "<RNA_DATA_PATH>",
        'riboseq_samples': riboseq_samples,
        'rnaseq_samples': rnaseq_samples,
        'riboseq_biological_replicates': ribo_condition_groups,
        'rnaseq_biological_replicates': rna_condition_groups
    }

    with open(args.out, 'w') as out:
        yaml.dump(config, out, default_flow_style=False)

if __name__ == '__main__':
    main()
