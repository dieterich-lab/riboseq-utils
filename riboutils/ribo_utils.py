import logging
import os
import pandas as pd

import riboutils.ribo_filenames as filenames

logger = logging.getLogger(__name__)


class _return_key_dict(dict):
    def __missing__(self,key):
        return key


###
#   The following labels are used to group similar ORF types.
###
orf_type_labels_mapping = {
    'canonical': ['canonical'],
    'canonical_variant': ['canonical_extended', 'canonical_truncated'],
    'five_prime': ['five_prime'],
    'three_prime': ['three_prime'],
    'noncoding': ['noncoding'],
    'novel': ['novel'],
    'other': [
        'five_prime_overlap', 
        'suspect_overlap', 
        'three_prime_overlap', 
        'within'
    ],
    'novel_overlap': [
        'novel_canonical', 
        'novel_canonical_extended', 
        'novel_canonical_truncated', 
        'novel_five_prime', 
        'novel_three_prime', 
        'novel_noncoding', 
        'novel_novel', 
        'novel_five_prime_overlap',
        'novel_suspect_overlap', 
        'novel_three_prime_overlap', 
        'novel_within'
    ]
}

orf_type_labels_reverse_mapping = {
    v:k for k, l in orf_type_labels_mapping.items() for v in l
}

orf_type_labels_display_name_map = {
    'canonical': "Canonical",
    'canonical_variant': "Canonical variant",
    'five_prime': "uORF",
    'three_prime': "dORF",
    'noncoding': "ncRNA",
    'novel': "de novo only",
    'other': "Other",
    'novel_overlap': "de novo overlap"
}

orf_type_display_name_map = {
    'canonical': "Canonical", 
    'canonical_extended': "Canonical extended", 
    'canonical_truncated': "Canonical truncated", 
    'five_prime': "uORF", 
    'three_prime': "dORF", 
    'noncoding': "ncRNA", 
    'five_prime_overlap': "uORF overlap",
    'suspect_overlap': "Suspect", 
    'three_prime_overlap': "dORF overlap", 
    'within': "Within",
    "novel": "de novo only",
    'novel_canonical_extended': "de novo canonical extended", 
    'novel_five_prime': "de novo uORF", 
    'novel_three_prime': "de novo dORF", 
    'novel_noncoding': "de novo ncRNA", 
    'novel_five_prime_overlap': "de novo uORF overlap",
    'novel_suspect_overlap': "de novo suspect", 
    'novel_three_prime_overlap': "de novo dORF overlap", 
    'novel_within': "de novo within",
}

orf_type_labels = list(orf_type_labels_mapping.keys())
orf_types = list(orf_type_display_name_map.keys())

###
#   The following functions are helpful for parsing information out of the identifiers.
###
def get_transcript_id(orf_id, sep="_"):
    return orf_id.split(sep)[0]

def get_all_transcript_ids(orfs, sep="_", num_cpus=1, progress_bar=False):
    import misc.parallel as parallel

    transcript_ids = parallel.apply_parallel_iter(  orfs['id'],
                                                    num_cpus,
                                                    get_transcript_id,
                                                    sep,
                                                    progress_bar=progress_bar)

    return transcript_ids

###
#   The following functions are all used for parsing replicates, etc., from the config file.
###

def get_sample_reverse_map(config):
    """ Extract a mapping from riboseq and rnaseq samples to conditions. """
    reverse_map = _return_key_dict()

    riboseq_reverse_map = get_riboseq_replicates_reverse_map(config)
    rnaseq_reverse_map = get_rnaseq_replicates_reverse_map(config)

    reverse_map.update(riboseq_reverse_map)
    reverse_map.update(rnaseq_reverse_map)
    return reverse_map

def get_riboseq_replicates(config):
    if 'riboseq_biological_replicates' in config:
        if config['riboseq_biological_replicates'] is not None:
            msg = "Found 'riboseq_biological_replicates' key in config file"
            logger.info(msg)
            return config['riboseq_biological_replicates']
        
    msg = ("Did not find 'riboseq_biological_replicates' key in config file. "
            "Using each 'riboseq_sample' as a single-condition replicate.")
    logger.info(msg)

    # create a dictionary mapping from the sample name to asingle-element list
    ret = {
        name: [name] for name, sample in config['riboseq_samples'].items()
    }
    return ret

def get_riboseq_replicates_reverse_map(config):
    """ Extract a mapping from sample to condition. """
    riboseq_replicates = get_riboseq_replicates(config)
    reverse_map = {
        v:k for k, l in riboseq_replicates.items() for v in l
    }
    return reverse_map


def get_field_condition_name_map(config):
    """ Extract a mapping from riboseq and rnaseq conditions to pretty names.
    """
    condition_name_map = _return_key_dict()

    riboseq_map = get_riboseq_condition_name_map(config)
    rnaseq_map = get_rnaseq_condition_name_map(config)

    condition_name_map.update(riboseq_map)
    condition_name_map.update(rnaseq_map)
    return condition_name_map


def get_riboseq_condition_name_map(config):
    """ Extract the pretty names for the riboseq replicates, if they are given
    in the config. All other names are returned unchanged.

    This is based on the 'riboseq_condition_name_map' key.
    """

    riboseq_condition_name_map = _return_key_dict()

    if 'riboseq_condition_name_map' in config:
        riboseq_condition_name_map.update(config['riboseq_condition_name_map'])

    return riboseq_condition_name_map


def get_rnaseq_condition_name_map(config):
    """ Extract the pretty names for the rnaseq conditions, if they are given
    in the config. All other names are returned unchanged.

    This is based on the 'rnaseq_condition_name_map' key.
    """

    rnaseq_condition_name_map = _return_key_dict()

    if 'rnaseq_condition_name_map' in config:
        rnaseq_condition_name_map.update(config['rnaseq_condition_name_map'])

    return rnaseq_condition_name_map



def get_rnaseq_replicates(config):
    if 'rnaseq_biological_replicates' in config:
        if config['rnaseq_biological_replicates'] is not None:
            msg = "Found 'rnaseq_biological_replicates' key in config file"
            logger.info(msg)

            return config['rnaseq_biological_replicates']
        
    msg = ("Did not find 'rnaseq_biological_replicates' key in config file. "
            "Using each 'rnaseq_sample' as a single-condition replicate.")
    logger.info(msg)
    
    # create a dictionary mapping from the sample name to asingle-element list
    ret = {
        name: [name] for name, sample in config['rnaseq_samples'].items()
    }
    return ret


def get_rnaseq_replicates_reverse_map(config):
    """ Extract a mapping from sample to condition. """
    rnaseq_replicates = get_rnaseq_replicates(config)
    reverse_map = {
        v:k for k, l in rnaseq_replicates.items() for v in l
    }
    return reverse_map


def get_matching_conditions(config):
    if 'matching_conditions' in config:
        if config['matching_conditions'] is not None:
            msg = "Found 'matching_conditions' key in config file"
            logger.info(msg)

            return config['matching_conditions']
        
    msg = ("Did not find 'matching_conditions' key in config file. Using "
            "riboseq and rnaseq conditions (biological_replicate entries) "
            "as matching conditions.")
    logger.info(msg)
    
    # otherwise, get the replicates and match key names
    riboseq_replicates = get_riboseq_replicates(config)
    rnaseq_replicates = get_rnaseq_replicates(config)
    
    matching_conditions = {
        x: [x, x] for x in riboseq_replicates if x in rnaseq_replicates
    }
    
    return matching_conditions

def get_matching_condition_and_replicates(condition:str, config:dict, 
        names_only:bool=False, raise_on_error:bool=True):
    """ Retrieve the matching ribo and rnaseq conditions for the given
    matching condition name from config.

    Parameters
    ----------
    condition: string
        the name of the "matching" condition

    config: dict
        the configuration dictionary

    names_only: bool
        whether to return only the matching ribo and rnaseq conditions

    raise_on_error: bool
        whether to raise an error or issue a warning message when values are 
        misssing

    Returns
    -------
    None:
        if raise_on_error is False and any keys are not found

    ... otherwise ...
    ribo_condition, rna_condition: strings
        the name of the respective conditions for this "matching" condition

    ribo_replicates, rna_replicates: list of strings
        the replicates for the respective conditions
    """
    # make sure the matching_condition exists
    matching_conditions = get_matching_conditions(config)
    if condition not in matching_conditions:
        msg = ("[ribo_utils.get_matching_condition_and_replicates]: Could not "
            "find \"{}\" in matching_conditions. Please ensure the name is "
            "spelled correctly.".format(condition))

        if raise_on_error:
            raise ValueError(msg)
        else:
            logger.warning(msg)
            return None

    # also, make sure the condition is in both of the replicate lists
    cond = matching_conditions[condition]

    if len(cond) != 2:
        msg = ("[ribo_utils.get_matching_condition_and_replicates]: A set of "
            "matching conditions is ill-formed. Each set of matching "
            "conditions must be a 2-tuple. This first condition should be the "
            "riboseq condition, and the second should be the rnaseq "
            "condition. '{}: {}'".format(condition, cond))

        if raise_on_error:
            raise ValueError(msg)
        else:
            logger.warning(msg)
            return None

    ribo_condition = cond[0]
    rna_condition = cond[1]

    riboseq_biological_replicates = get_riboseq_replicates(config)
    rnaseq_biological_replicates = get_rnaseq_replicates(config)

    if ribo_condition not in riboseq_biological_replicates:
        msg = ("[ribo_utils.get_matching_condition_and_replicates]: The "
            "riboseq condition '{}' is not present in the "
            "'riboseq_biological_replicates'.".format(ribo_condition))

        if raise_on_error:
            raise ValueError(msg)
        else:
            logger.warning(msg)
            return None

    if rna_condition not in rnaseq_biological_replicates:
        msg = ("[ribo_utils.get_matching_condition_and_replicates]: The rna "
            "condition '{}' is not present in the "
            "'rnaseq_biological_replicates'.".format(rna_condition))
        
        if raise_on_error:
            raise ValueError(msg)
        else:
            logger.warning(msg)
            return None

    if names_only:
        return ribo_conditions, rna_conditions

    ribo_replicates = riboseq_biological_replicates[ribo_condition]
    rna_replicates = rnaseq_biological_replicates[rna_condition]

    return ribo_condition, rna_condition, ribo_replicates, rna_replicates


def get_criterion_condition(condition, criterion, config):
    matching_conditions = get_matching_conditions(config)
    if condition not in matching_conditions:
        msg = ("[ribo_utils.get_criterion_condition]: Could not find '{}' in "
            "'matching_conditions".format(condition))
        raise ValueError(msg)

    ribo, rna = matching_conditions[condition]

    if criterion == "ribo":
        return ribo
    elif criterion == "rna":
        return rna
    else:
        msg = ("[ribo_utils.get_criterion_condition]: The criterion '{}' is "
            "not a valid criterion.".format(criterion))
        raise ValueError(msg)

def get_riboseq_cell_type_samples(config):
    if 'riboseq_cell_type_samples' in config:
        if config['riboseq_cell_type_samples'] is not None:
            msg = "Found 'riboseq_cell_type_samples' key in config file"
            logger.info(msg)
            return config['riboseq_cell_type_samples']

    msg = ("Did not find 'riboseq_cell_type_samples' key in config file. Using "
            "riboseq conditions (biological_replicate entries) as the cell types")
    logger.info(msg)

    riboseq_replicates = get_riboseq_replicates(config)
    cell_type_samples = {
        x: [x] for x in riboseq_replicates
    }
    return cell_type_samples


def get_rnaseq_cell_type_samples(config):
    if 'rnaseq_cell_type_samples' in config:
        if config['rnaseq_cell_type_samples'] is not None:
            msg = "Found 'rnaseq_cell_type_samples' key in config file"
            logger.info(msg)
            return config['rnaseq_cell_type_samples']

    msg = ("Did not find 'rnaseq_cell_type_samples' key in config file. Using "
            "riboseq conditions (biological_replicate entries) as the cell types")
    logger.info(msg)

    rnaseq_replicates = get_rnaseq_replicates(config)
    cell_type_samples = {
        x: [x] for x in rnaseq_replicates
    }
    return cell_type_samples

def get_sample_name_map(config):
    """ Extract the mapping from the '{ribo,rna}seq_sample_name_map', or create
    a default one for all samples without an entry.
    """

    sample_name_map = _return_key_dict()

    if 'riboseq_sample_name_map' in config:
        sample_name_map.update(config['riboseq_sample_name_map'])

    if 'rnaseq_sample_name_map' in config:
        sample_name_map.update(config['rnaseq_sample_name_map'])

    return sample_name_map

def get_condition_name_map(config):
    """ Extract the mapping from the 'condition_name_map' and create a default
    one for all conditions without an entry.
    """
    condition_name_map = _return_key_dict()
    if 'condition_name_map' in config:
        condition_name_map.update(config['condition_name_map'])

    return condition_name_map

def filter_condition_pairs(config, allowed_conditions):
    """ Create an iterator which yields only condition pairs for which both 
    conditions appear in the allowed_conditions.

    Parameters
    ----------
    config: dict (presumably loaded from a yaml config file)
        A configuration dictionary which *must* include comparison_conditions

    allowed_conditions: sequence or None
        The conditions we care about. If None or the length is 0, then none of
        the condition pairs are filtered (all are "yield"ed).

    Yields
    ------
    condition_pair: 2-tuple of strings
        The next condition pair which meets the filtering criteria
    """
    import misc.utils as utils

    condition_pairs = config['comparison_conditions']

    if (allowed_conditions is not None) and (len(allowed_conditions) > 0):
        allowed_conditions = set(allowed_conditions)
    else:
        allowed_conditions = set(utils.flatten_lists(condition_pairs))

    for cp in condition_pairs:
        if (cp[0] in allowed_conditions) and (cp[1] in allowed_conditions):
            yield cp



###
#
# This function is used to extract the lengths and offsets which count as
# "periodic," based on the values in the config file.
#
###

default_min_metagene_profile_count = 1000
default_min_metagene_bf_mean = 5
default_max_metagene_bf_var = None
default_min_metagene_bf_likelihood = 0.5

def get_periodic_lengths_and_offsets(config, name, do_not_call=False,
        isoform_strategy=None, is_unique=True):

    """ This function applies a set of filters to metagene profiles to select those
        which are "periodic" based on the read counts and Bayes factor estimates.

        First, the function checks if the configuration file sets the 
        'use_fixed_lengths' flag is set. If so, then the specified lengths and
        offsets are returned.

        Otherwise, the function opens the appropriate file and extracts the filter
        values from the configuration file. In particular, it looks for the
        following keys:

        min_metagene_profile_count (float) : the minimum number of reads for a 
            particular length in the filtered genome profile. Read lengths with 
            fewer than this number of reads will not be used. default: 1000

        min_metagene_bf_mean (float) : if max_metagene_profile_bayes_factor_var 
            is not None, then this is taken as a hard threshold on the estimated 
            Bayes factor mean. If min_metagene_profile_bayes_factor_likelihood is 
            given, then this is taken as the boundary value; that is, a profile is
            "periodic" if:

                    [P(bf > min_metagene_bf_mean)] > min_metagene_bf_likelihood

            If both max_metagene_bf_var and min_metagene_bf_likelihood are None, 
            then this is taken as a hard threshold on the mean for selecting 
            periodic read lengths.

            If both max_metagene_bf_var and min_metagene_bf_likelihood are given, 
            then both filters will be applied and the result will be the intersection.

        max_metagene_bf_var (float) : if given, then this is taken as a hard threshold
            on the estimated Bayes factor variance. default: None (i.e., this filter
            is not used)

        min_metagene_bf_likelihood (float) : if given, then this is taken a threshold
            on the likelihood of periodicity (see min_metagene_bf_mean description
            for more details). default: 0.5

        Parameters
        ----------
        config: dictionary
            the configuration information(see description)

        name: string
            the name of the dataset in question

        do_not_call: bool
            whether the metagene bf file should exist. If false, then dummy
            values are returned (and a warning message is printed).

        isoform_strategy: string
            which strategy is used to select isoforms (relevant for B-tea only)

        is_unique: bool
            whether only unique reads are used in the files

        Returns
        -------
        lengths: list of strings
            all of the periodic read lengths

        offsets: list of strings
            the corresponding P-site offsets for the read lengths
    """
    import numpy as np
    import scipy.stats

    # check if we specified to just use a fixed offset and length
    if 'use_fixed_lengths' in config:
        lengths = [str(l) for l in config['lengths']]
        offsets = [str(o) for o in config['offsets']]

        return (lengths, offsets)

    # filter out the lengths which do not satisfy the quality thresholds
    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", default_min_metagene_profile_count)

    min_bf_mean = config.get(
        "min_metagene_bf_mean", default_min_metagene_bf_mean)

    max_bf_var = config.get(
        "max_metagene_bf_var", default_max_metagene_bf_var)
        
    min_bf_likelihood = config.get(
        "min_metagene_bf_likelihood", default_min_metagene_bf_likelihood)

    note_str = config.get('note', None)

    periodic_offsets = filenames.get_periodic_offsets(
        config['riboseq_data'],
        name, 
        is_unique=is_unique,
        isoform_strategy=isoform_strategy,
        note=note_str
    )
    
    if not os.path.exists(periodic_offsets):
        msg = ("The periodic offsets file does not exist. Please ensure the "
            "select-periodic-offsets script completed successfully or specify "
            "the \"use_fixed_lengths\", \"lengths\", and \"offsets\" values "
            "in the configuration file. '{}'".format(periodic_offsets))

        if do_not_call:
            msg = msg +  ("\nThe --do-not-call flag was given, so a \"dummy\" "
                "default length (29) and offset (12) will be used to check "
                "the remaining calls.\n")

            logger.warning(msg)

            offsets = ["12"]
            lengths = ["29"]
            return (lengths, offsets)
        else:
            raise FileNotFoundError(msg)
    
    offsets_df = pd.read_csv(periodic_offsets)


    # we always use the count filter
    m_count = offsets_df['highest_peak_profile_sum'] > min_metagene_profile_count

    # which bf mean/variance filters do we use? 
    m_bf_mean = True
    m_bf_var = True
    m_bf_likelihood = True

    if max_bf_var is not None:
        m_bf_mean = offsets_df['highest_peak_bf_mean'] > min_bf_mean
        m_bf_var = offsets_df['highest_peak_bf_var'] < max_bf_var

        msg = ("Using the mean and variance filter. min_mean: {}, max_var: {}"
            .format(min_bf_mean, max_bf_var))
        logger.debug(msg)

    if min_bf_likelihood is not None:
        # first, calculate the likelihood that the true BF is greater than m_bf_mean

        # the likelihood that BF>min_mean is 1-cdf(estimated_mean, estimated_var)

        # scipy parameterizes the normal using the std, so use sqrt(var)

        likelihood = 1-scipy.stats.norm.cdf(min_bf_mean, offsets_df['highest_peak_bf_mean'], 
            np.sqrt(offsets_df['highest_peak_bf_var']))

        nans = np.isnan(likelihood)
        num_nans = sum(nans)
        num_predictions = len(likelihood)

        max_likelihood = max(likelihood[~nans])

        # now filter
        m_bf_likelihood = likelihood > min_bf_likelihood

    if (max_bf_var is None) and (min_bf_likelihood is None):
        m_bf_mean = bf['highest_peak_bf_mean'] > min_bf_mean

    filtered_periodic_offsets = offsets_df[m_count & m_bf_mean & m_bf_var & m_bf_likelihood]

    offsets = filtered_periodic_offsets['highest_peak_offset']
    lengths = filtered_periodic_offsets['length']

     
    if len(lengths) == 0:
        msg = ("The periodic offsets file was found, but no periodic lengths "
            "were found. Please ensure the select-periodic-offsets script "
            "completed successfully or specify the \"use_fixed_lengths\", "
            "\"lengths\", and \"offsets\" values in the configuration file. "
            "'{}'".format(periodic_offsets))

        if do_not_call:
            msg = msg +  ("\nThe --do-not-call flag was given, so a \"dummy\" "
                "default length (29) and offset (12) will be used to check "
                "the remaining calls.\n")

            logger.warning(msg)

            offsets = ["12"]
            lengths = ["29"]
            return (lengths, offsets)
        else:
            raise ValueError(msg)


    # offsets must be positive
    offsets = [str(-1*int(o)) for o in offsets]
    lengths = [str(int(l)) for l in lengths]

    return (lengths, offsets)

###
#   This function extracts the p-sites from alignments, given the offsets
#   and periodic read lengths.
###
def get_p_sites(bam_file, periodic_lengths, offsets):
    """ Given a bam file of mapped riboseq reads, this function filters
        out the reads of non-periodic length, adjusts the start and end
        positions based on strand, and then shifts the remaining reads
        based on the length-specific offset.
        
        Args:
            bam_file (string) : the path to the mapped riboseq reads
            
            periodic_lengths (list-like) : a list of lengths to keep
            
            offsets (list-like) : the distance to shift each read of the
                respective length. the order here must match that in
                periodic_lengths
                
        Returns:
            pd.DataFrame : a data frame containing the transformed reads,
                sorted by chrom and start

        Imports:
            sys
            numpy
            pandas
            tqdm
            pysam
            misc.bio
    """
    import sys
    import numpy as np
    import pandas as pd
    import tqdm

    import pysam
    import misc.bio as bio

    msg = "Reading BAM file"
    logger.info(msg)

    bam = pysam.AlignmentFile(bam_file)
    alignments = bam.fetch()
    num_alignments = bam.count()

    logger.info("Processing alignments")

    lengths = np.zeros(num_alignments, dtype=int)
    starts = np.zeros(num_alignments, dtype=int)
    ends = np.zeros(num_alignments, dtype=int)
    seqs = [""] * num_alignments
    strands = ["+"] * num_alignments
    fractions = np.zeros(num_alignments, dtype=float)

    al_iter = tqdm.tqdm(alignments, leave=True, file=sys.stdout, total=num_alignments)
    for i, a in enumerate(al_iter):
        starts[i] = a.reference_start
        ends[i] = a.reference_end
        lengths[i] = a.qlen
        seqs[i] = a.reference_name

        if a.is_reverse:
            strands[i] = "-"

    # The data frame will later be converted to BED6, so put the fields in the
    # correct order.
    map_df = pd.DataFrame()
    map_df['seqname'] = seqs
    map_df['start'] = starts
    map_df['end'] = ends
    map_df['id'] = "."
    map_df['score'] = "."
    map_df['strand'] = strands
    map_df['length'] = lengths

    msg = "Filtering reads by length"
    logger.info(msg)
    
    # now, filter based on lengths
    m_length = map_df['length'].isin(periodic_lengths)
    map_df = map_df[m_length]

    # now, we need to update the starts and ends based on the strand
    msg = "Updating coordinates based on offsets"
    logger.info(msg)

    # if the strand is positive, the end is start+1
    # if the strand is negative, the start is end-1
    m_positive = map_df['strand'] == '+'
    m_negative = map_df['strand'] == '-'
    
    # first, shift in the appropriate direction
    for i in range(len(periodic_lengths)):
        length = periodic_lengths[i]
        offset = offsets[i]
        
        m_length = map_df['length'] == length
        
        # adjust the start of forward strand
        map_df.loc[m_positive & m_length, 'start'] = (
                map_df.loc[m_positive & m_length, 'start'] + offset)
        
        # adjust the ends of negative strand
        map_df.loc[m_negative & m_length, 'end'] = (
                map_df.loc[m_negative & m_length, 'end'] - offset)

    # finally, we only care about the 5' end of the read, so discard everything else
    msg = "Discarding 3' end of reads"
    logger.info(msg)
    
    map_df.loc[m_positive, 'end'] = map_df.loc[m_positive, 'start'] + 1
    map_df.loc[m_negative, 'start'] = map_df.loc[m_negative, 'end'] - 1

    # now sort everything
    msg = "Sorting reads by coordinates"
    logger.info(msg)
    
    map_df = map_df.sort_values(['seqname', 'start'])

    # and we only want the BED6 fields
    map_df = map_df[bio.bed6_field_names]
    
    return map_df

###
#   This function smoothes the profiles, frame-by-frame
###

default_fraction = 0.2
default_reweighting_iterations = 0

def smooth_profile(profile, reweighting_iterations=default_reweighting_iterations, 
        fraction=default_fraction):

    """ This function smoothes the given ORF profile using the frame-specific
        approach. It assumes the profile is a dense numpy array and that any
        filtering due to differences of counts in reading frames, lengths, etc.,
        has already been performed.

        Please see the statsmodels.api.nonparametric.lowess documentation for
        more information about reweighting_iterations and fraction.

        Args:
            profile (np.array of numbers): an array containing the observed
                ORF profile. In principle, this could already be normalized.

            reweighting_iterations (int): the number of reweighting iterations

            fraction (float): the percentage of the signal to use for smooothing

        Returns:
            np.array: the smoothed profile

        Imports:
            statsmodels.api.nonparametric.lowess
    """
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess
    import numpy as np


    smoothed_profile = np.zeros_like(profile)

    # split the signal based on frame
    x_1 = profile[0::3]
    x_2 = profile[1::3]
    x_3 = profile[2::3]
    exog = np.arange(len(x_1))

    # x_1
    endog = x_1
    smoothed_x_1 = lowess(endog, exog, is_sorted=True, return_sorted=False, 
        it=reweighting_iterations, frac=fraction)
    
    # x_2
    endog = x_2
    smoothed_x_2 = lowess(endog, exog, is_sorted=True, return_sorted=False, 
        it=reweighting_iterations, frac=fraction)
    
    # x_3
    endog = x_3
    smoothed_x_3 = lowess(endog, exog, is_sorted=True, return_sorted=False, 
        it=reweighting_iterations, frac=fraction)
    
    smoothed_profile[0::3] = smoothed_x_1
    smoothed_profile[1::3] = smoothed_x_2
    smoothed_profile[2::3] = smoothed_x_3

    return smoothed_profile


###
#
# This function extracts all ORFs which count as "translated", according to the
# values in the config file.
#
###

default_min_profile = None
default_min_bf_mean = 5
default_max_bf_var = None
default_min_bf_likelihood = None
default_min_length = 0
default_chisq_alpha = 0.01

def get_base_filter(bf, min_profile=default_min_profile, min_length=default_min_length):
    """ This function extracts the ORFs from the BF dataframe which meet the
        minimum requirements to be considered for prediction. Namely, these
        requirements are:
        
            * The minimum sum across all reading frames exceeds the specified minimum
            * The length exceeds the specified minimum length
            * The number of reads in the first reading frame exceeds the number in
                either of the other two reading frames (though not necessarily the
                other two reading frames combined)

        Args:
            bf (pd.DataFrame): a data frame containing the relevant ORF information

            min_signal (int) : the minimum sum across all reading frames to consider
                an ORF as translated
            
            min_length (int) : the minimum length of ORF to consider

        Returns:
            boolean mask: a mask of the input data frame indicating all ORFs which
                meet the filtering criteria
    """
    
    if min_profile is None:
        m_profile = bf['profile_sum'] > 0
    else:
        m_profile = bf['profile_sum'] > min_profile

    m_length = bf['orf_len'] > min_length
    m_x1_gt_x2 = bf['x_1_sum'] > bf['x_2_sum']
    m_x1_gt_x3 = bf['x_1_sum'] > bf['x_3_sum']

    m_base = m_profile & m_length & m_x1_gt_x2 & m_x1_gt_x3
    return m_base

def get_bf_filter(bf, min_bf_mean=default_min_bf_mean, 
                    max_bf_var=default_max_bf_var,
                    min_bf_likelihood=default_min_bf_likelihood):

    """ This function applies filters to the Bayes factor estimates to find all
        ORFs which should be predicted as translated. This does not consider the
        length and profile sums, so this filter would need to be combined with
        the get_base_filter filter to find the true set of predicted ORFs.

        Args:
            bf (pd.DataFrame) : a data frame containing the relevant ORF information

            min_bf_mean (float) : if max_bf_var is not None, then this is taken
                as a hard threshold on the estimated Bayes factor mean. If
                min_bf_likelihood is given, then this is taken as the boundary
                value; that is, an ORF is "translated" if:

                    [P(bf > min_bf_mean)] > min_bf_likelihood

                If both max_bf_var and min_bf_likelihood are None, then this is
                taken as a hard threshold on the mean for selecting translated ORFs.

                If both max_bf_var and min_bf_likelihood are given, then both
                filters will be applied and the result will be the intersection.

            max_bf_var (float) : if given, then this is taken as a hard threshold
                on the estimated Bayes factor variance
            
            min_bf_likelihood (float) : if given, then this is taken a threshold
                on the likelihood of translation (see min_bf_mean description
                for more details)
        
        Returns:
            boolean mask: a mask of the input data frame indicating all ORFs which
                meet the filtering criteria

        Imports:
            numpy
            scipy.stats
    """
    import numpy as np
    import scipy.stats

    # which bf mean/variance filters do we use? 
    m_bf_mean = True
    m_bf_var = True
    m_bf_likelihood = True

    if max_bf_var is not None:
        m_bf_mean = bf['bayes_factor_mean'] > min_bf_mean
        m_bf_var = bf['bayes_factor_var'] < max_bf_var
    if min_bf_likelihood is not None:
        # first, calculate the likelihood that the true BF is greater than m_bf_mean

        # the likelihood that BF>min_mean is 1-cdf(estimated_mean, estimated_var)

        # scipy parameterizes the normal using the std, so use sqrt(var)

        loc = bf['bayes_factor_mean']
        scale = np.sqrt(bf['bayes_factor_var'])
        likelihood = 1-scipy.stats.norm.cdf(min_bf_mean, loc, scale)

        nans = np.isnan(likelihood)
        num_nans = sum(nans)
        num_predictions = len(likelihood)

        msg = "Num nans: {}, num predictions: {}".format(num_nans, num_predictions)
        logger.debug(msg)

        if num_nans != num_predictions:
            max_likelihood = max(likelihood[~nans])
            msg = "Maximum likelihood: {}".format(max_likelihood)
            logger.debug(msg)

        # now filter
        m_bf_likelihood = likelihood > min_bf_likelihood

    if (max_bf_var is None) and (min_bf_likelihood is None):
        m_bf_mean = bf['bayes_factor_mean'] > min_bf_mean

    return m_bf_mean & m_bf_var & m_bf_likelihood



def get_predicted_orfs(bf, min_signal=default_min_profile, 
                            min_length=default_min_length,
                            min_bf_mean=default_min_bf_mean, 
                            max_bf_var=default_max_bf_var,
                            min_bf_likelihood=default_min_bf_likelihood,
                            chisq_alpha=default_chisq_alpha,
                            select_longest_by_stop=True):
    """ This function applies a set of filters to ORFs to select those which
        are predicted as "translated." This function selects translated ORFs
        based on the Bayes factor estimates and the chi-square p-values. ORFs
        must pass all of the relevant features to be selected as "translated."
        Optionally, among all ORFs which share a stop codon, only the longest
        "translated" ORF is selected.

        Furthermore, for both BF and chi-square predictions, only ORFs which
        have more reads in the first reading frame than either of the other two
        will be selected as translated. (This is called the 'frame filter'
        below.)

        Args:
            bf (pd.DataFrame) : a data frame containing the relevant ORF information

            min_signal (int) : the minimum sum across all reading frames to consider
                an ORF as translated
            
            min_length (int) : the minimum length of ORF to consider

            min_bf_mean (float) : if max_bf_var is not None, then this is taken
                as a hard threshold on the estimated Bayes factor mean. If
                min_bf_likelihood is given, then this is taken as the boundary
                value; that is, an ORF is "translated" if:

                    [P(bf > min_bf_mean)] > min_bf_likelihood

                If both max_bf_var and min_bf_likelihood are None, then this is
                taken as a hard threshold on the mean for selecting translated ORFs.

                If both max_bf_var and min_bf_likelihood are given, then both
                filters will be applied and the result will be the intersection.

            max_bf_var (float) : if given, then this is taken as a hard threshold
                on the estimated Bayes factor variance

            min_bf_likelihood (float) : if given, then this is taken a threshold
                on the likelihood of translation (see min_bf_mean description
                for more details)

            chisq_alpha (float) : the significance value for selecting translated
                ORFs according to the chi-square test. This value is 
                Bonferroni-corrected based on the number of ORFs which meet the
                length, profile and frame filters.

            select_longest_by_stop (bool): if True, then the selected ORFs will
                be merged based on stop codons: only the longest translated ORF
                at each stop codon will be returned. Otherwise, all ORFs will
                be returned.

        Returns:
            all_orfs (pd.DataFrame) : all (longest) ORFs which meet the profile,
                 length, frame filters

            bf_longest_orfs (pd.DataFrame) : all longest ORFs which meet the
                profile, length, frame (min_bf_mean, max_bf_var, min_bf_likelihood) filters

            chisq_longest_orfs (pd.DataFrame) : all longest ORFs which meet the
                profile, length, frame, chisq_alpha filters

        Imports:
            misc.bio
            numpy
            scipy.stats

    """
    import misc.bio as bio
    import misc.bio_utils.bed_utils as bed_utils
    import numpy as np
    import scipy.stats

    msg = "Finding all ORFs with signal"
    logger.info(msg)

    m_base = get_base_filter(bf, min_signal, min_length)
    all_orfs = bf[m_base]
    
    
    # create the selected ORFs based on Bayes factor
    m_bf =  get_bf_filter(bf, min_bf_mean, max_bf_var, min_bf_likelihood)
    m_bf_predicted = m_base & m_bf
    bf_predicted_orfs = bf[m_bf_predicted]

    M = len(all_orfs)
    # for the bonferroni correction, we only correct for the number of tests we actually consider
    # that is, we only correct for orfs which pass the base filter
    corrected_significance_level = chisq_alpha / M

    msg = "Corrected significance level: {}".format(corrected_significance_level)
    logger.debug(msg)
    
    m_chisq_pval = bf['chi_square_p'] < corrected_significance_level
    m_chisq_predicted = m_base & m_chisq_pval
    chisq_predicted_orfs = bf[m_chisq_predicted]
    
    if select_longest_by_stop:
        all_orfs = bed_utils.get_longest_features_by_end(all_orfs)
        bf_predicted_orfs = bed_utils.get_longest_features_by_end(bf_predicted_orfs)
        chisq_predicted_orfs = bed_utils.get_longest_features_by_end(chisq_predicted_orfs)
    
    return (all_orfs, bf_predicted_orfs, chisq_predicted_orfs)
   
###
#   Defaults for b-tea scripts
###
default_perm_test_min_rpkm_mean = 1
default_perm_test_max_rpkm_var_power = 1

###
#   Field names for b-tea files
###

field_map = {
    "ribo": "Riboseq",
    "rna": "RNA-seq",
    "te": "Translational Efficiency"
}

fields = sorted(field_map.keys())
field_name_order = [field_map[f] for f in fields]

def get_field_name(field):
    """ This function maps from the field to a human-readable name.
    """
    
    return field_map[field]

mean_field_map = {
    "ribo": "ribo_abundance_mean_loc",
    "rna": "rna_abundance_mean_loc",
    "te": "log_translational_efficiency_loc"
}

var_field_map = {
 "ribo": "ribo_abundance_var_loc",
 "rna": "rna_abundance_var_loc",
 "te": "log_translational_efficiency_scale"
}

###
# The following functions are all related. They are used to estimate p-values
# for the KL-divergence values calculated for translational efficiency (only).
###
def get_basic_filter(kl, condition_1, condition_2, field):
    """ Mask kl to filter on the conditions and field. """
    m_condition_1 = kl['condition_1'] == condition_1
    m_condition_2 = kl['condition_2'] == condition_2
    m_field = kl['field'] == field
    m_basic = m_condition_1 & m_condition_2 & m_field
    return m_basic

def get_rpkm_mean_filter(kl, min_rpkm_mean):
    """ Mask kl to filter on the estimated means. """
    
    m_min_rpkm_mean_1 = kl['mean_1'] > min_rpkm_mean
    m_min_rpkm_mean_2 = kl['mean_2'] > min_rpkm_mean
    m_min_rpkm_mean = m_min_rpkm_mean_1 & m_min_rpkm_mean_2
    return m_min_rpkm_mean

def get_rpkm_var_power_filter(kl, max_rpkm_var_power):
    """ Mask kl to filter on the variances as a power of the means. """
    import numpy as np
    
    m_max_rpkm_var_1 = kl['var_1'] < np.power(kl['mean_1'], max_rpkm_var_power)
    m_max_rpkm_var_2 = kl['var_2'] < np.power(kl['mean_2'], max_rpkm_var_power)
    m_max_rpkm_var = m_max_rpkm_var_1 & m_max_rpkm_var_2
    return m_max_rpkm_var

def get_basic_and_rpkm_filter(kl, condition_1, condition_2, field, 
        min_rpkm_mean, max_rpkm_var_power):
    """ Mask kl using all of the indicated filters. This handles TE as the 
    combination of both riboseq and rnaseq.
    """

    if field == "te":
        # first, get the genes which meet the rpkm requirements
        m_ribo = get_basic_and_rpkm_filter(
            kl, 
            condition_1, 
            condition_2, 
            "ribo", 
            min_rpkm_mean, 
            max_rpkm_var_power
        )
        
        m_rna = get_basic_and_rpkm_filter(
            kl, 
            condition_1, 
            condition_2, 
            "rna", 
            min_rpkm_mean, 
            max_rpkm_var_power
        )

        # find the gene ids that meet both filters
        ribo_gene_ids = set(kl.loc[m_ribo, 'gene_id'].unique())
        rna_gene_ids = set(kl.loc[m_rna, 'gene_id'].unique())
        gene_ids = ribo_gene_ids & rna_gene_ids

        # get all te rows for these conditions
        m_basic = get_basic_filter(kl, condition_1, condition_2, "te")

        # and only keep the genes which met both rpkm requirements
        m_gene_ids = kl['gene_id'].isin(gene_ids)
        m_all = m_basic & m_gene_ids

    else:
        m_basic = get_basic_filter(kl, condition_1, condition_2, field)
        m_min_rpkm_mean = get_rpkm_mean_filter(kl, min_rpkm_mean)
        m_max_rpkm_var = get_rpkm_var_power_filter(kl, max_rpkm_var_power)
        m_all = m_basic & m_min_rpkm_mean & m_max_rpkm_var
    return m_all

###
# These functions are all based on the old "wide" data frame format. Thus, they
# have all been deprecated.
###

mean_format_map = {
    "te": "log_translational_efficiency_loc_{1}",
    "ribo": "{}_abundance_mean_loc_{}",
    "rna": "{}_abundance_mean_loc_{}"
}

var_format_map = {
    "te": "log_translational_efficiency_scale_{1}",
    "ribo": "{}_abundance_var_loc_{}",
    "rna": "{}_abundance_var_loc_{}"
}

# decorator to raise the deprecated warning
def ribo_deprecated(func):
    """ Issue a warning that the given function uses the "wide" df format and
        should be replaced with the easier to work with "long" format.
    """
    def wrapper(*args, **kwargs):
        msg = ("[ribo_utils.{}]: This function has been deprecated. It uses "
            "the old \"wide\" df format. Please replace it with the "
            "respective \"long\" df format function.".format(func.__name__))
        logger.warning(msg)

        return func(*args, **kwargs)
    return wrapper

@ribo_deprecated
def get_mean_var_column_names(field, condition):
    """ This function returns the name of the columns containing the mean and
        variance for the given field and condition.

        Parameters
        ----------
        field : string
            The name of the field in question. Valid values are:
                * te
                * ribo
                * rna

        condition : string
            The name of the condition (e.g., "sham.cm")

        Returns
        -------
        mean_column : string
            The name of the column containing the means for this field

        var_column : string
            The name of the column containing the variances for this field
    """
    mean_field = mean_format_map[field].format(field, condition)
    var_field = var_format_map[field].format(field, condition)

    return (mean_field, var_field)

kl_format_map = {
    "te": "log_translational_efficiency_{}_{}_kl_divergence",
    "ribo": "ribo_abundance_{}_{}_kl_divergence",
    "rna": "rna_abundance_{}_{}_kl_divergence"
}

pvalue_format_map = {
    "te": "log_translational_efficiency_{}_{}_pvalue",
    "ribo": "ribo_abundance_{}_{}_pvalue",
    "rna": "rna_abundance_{}_{}_pvalue"
}

@ribo_deprecated
def get_kl_pvalue_column_name(field, condition_1, condition_2):
    """ This function returns the names of the columns containing the estimated
        KL-divergence and pvalues for the two conditions and field.

    Parameters
    ----------
    field : string
        The name of the field in question. Valid values are:
            * te
            * ribo
            * rna

    condition_{1,2} : string
        The name of the condition (e.g., "sham.cm")

    Returns
    -------
    kl_column : string
        The name of the column containing the KL-divergence for this field

    pvalue_column : string
        The name of the column containing the means-values for this field
    """
    kl_field = kl_format_map[field].format(condition_1, condition_2)
    pvalue_field = pvalue_format_map[field].format(condition_1, condition_2)

    return (kl_field, pvalue_field)

significant_pvalue_format_map = {
    "te": "significant_te_{}_{}",
    "ribo": "significant_ribo_{}_{}",
    "rna": "significant_rna_{}_{}"
}

@ribo_deprecated
def get_significant_pvalue_column_name(field, condition_1, condition_2):
    """ Column name indicating the specified estimates significantly differ

    Parameters
    ----------
    field : string
        The name of the field in question. Valid values are:
            * te
            * ribo
            * rna

    condition_{1,2} : string
        The name of the conditions (e.g., "sham.cm")

    Returns
    -------
    significant_column: string
        Name of the column indicating significance
    """
    sig_pval_col = significant_pvalue_format_map[field]
    sig_pval_col = sig_pval_col.format(condition_1, condition_2)
    return sig_pval_col

@ribo_deprecated
def get_micropeptide_overlap_column_name(condition):
    """ Column name indicating an overlap with a micropeptide

    Parameters
    ----------
    condition: string
        The name of the condition (e.g., "sham.cm")

    Returns
    -------
        Name of the column indicating an overlap
    """
    return "has_micropeptide_overlap_{}".format(condition)

log_fold_change_map = {
    "te": "log_translational_efficiency_{}_{}_log_fold_change",
    "ribo": "ribo_abundance_{}_{}_log_fold_change",
    "rna": "rna_abundance_{}_{}_log_fold_change"
}

@ribo_deprecated
def get_log_fold_change_field_name(field, condition_1, condition_2):
    lfc_field = log_fold_change_map[field].format(condition_1, condition_2)
    return lfc_field

@ribo_deprecated
def get_log_fold_changes(df, condition_pairs):
    """ This function creates a new data frame which includes all of the log
        fold changes (TE, riboseq and RNA-seq) for each of the condition
        pairs in the given list.

    The returned data frame could be joined to the original df with a
    command like:

    pd.concat([df, log_fold_changes_df], axis=1)

    Parameters
    ----------
    df : pd.DataFrame
        A data frame containing the "mean" fields

    condition_pairs : list of 2-tuple-likes of strings
        The pairs of conditions for which the log fold changes will be
        included in the returns data frame

    Returns
    -------
    log_fold_changes_df : pd.DataFrame
        A data frame containing all of the requested log fold changes
    """
    import numpy as np
    import pandas as pd
    
    log_fold_changes_df = pd.DataFrame()

    for (condition_1, condition_2) in condition_pairs:
        
        field = 'te'
        field_1 = mean_format_map[field].format(field, condition_1)
        field_2 = mean_format_map[field].format(field, condition_2)
        lfc_field = log_fold_change_map[field].format(condition_1, condition_2)
        log_fold_changes_df[lfc_field] = df[field_2] - df[field_1]
        
        field = 'ribo'
        field_1 = mean_format_map[field].format(field, condition_1)
        field_2 = mean_format_map[field].format(field, condition_2)
        lfc_field = log_fold_change_map[field].format(condition_1, condition_2)
        log_fold_changes_df[lfc_field] = np.log(df[field_2]) - np.log(df[field_1])
        
        field = 'rna'
        field_1 = mean_format_map[field].format(field, condition_1)
        field_2 = mean_format_map[field].format(field, condition_2)
        lfc_field = log_fold_change_map[field].format(condition_1, condition_2)
        log_fold_changes_df[lfc_field] = np.log(df[field_2]) - np.log(df[field_1])
               
    return log_fold_changes_df

@ribo_deprecated
def get_variance_power_filter(kl_df, condition_1, condition_2, field, power=0.5):
    import numpy as np

    # first, get the field names for which we want significances
    if field == "log_translational_efficiency":
        
        # filter by both rna_abundance and ribo_abundance in both samples
        ribo_var_1_f = "ribo_abundance_var_loc_{}".format(condition_1)
        ribo_var_2_f = "ribo_abundance_var_loc_{}".format(condition_2)

        rna_var_1_f = "rna_abundance_var_loc_{}".format(condition_1)
        rna_var_2_f = "rna_abundance_var_loc_{}".format(condition_2)

        # filter by both rna_abundance and ribo_abundance in both samples
        ribo_mean_1_f = "ribo_abundance_mean_loc_{}".format(condition_1)
        ribo_mean_2_f = "ribo_abundance_mean_loc_{}".format(condition_2)

        rna_mean_1_f = "rna_abundance_mean_loc_{}".format(condition_1)
        rna_mean_2_f = "rna_abundance_mean_loc_{}".format(condition_2)

        m_ribo_1 = abs(kl_df[ribo_var_1_f]) < np.power(abs(kl_df[ribo_mean_1_f]), power)
        m_ribo_2 = abs(kl_df[ribo_var_2_f]) < np.power(abs(kl_df[ribo_mean_2_f]), power)
        
        m_rna_1 = abs(kl_df[rna_var_1_f]) < np.power(abs(kl_df[rna_mean_1_f]), power)
        m_rna_2 = abs(kl_df[rna_var_2_f]) < np.power(abs(kl_df[rna_mean_2_f]), power)
        
        m_filter = (m_ribo_1 & m_ribo_2 & m_rna_1 & m_rna_2)

    else:
        var_1_f = "{}_var_loc_{}".format(field, condition_1)
        var_2_f = "{}_var_loc_{}".format(field, condition_2)

        mean_1_f = "{}_mean_loc_{}".format(field, condition_1)
        mean_2_f = "{}_mean_loc_{}".format(field, condition_2)


        # also get the filter
        m_1 = abs(kl_df[var_1_f]) < np.power(abs(kl_df[mean_1_f]), power)
        m_2 = abs(kl_df[var_2_f]) < np.power(abs(kl_df[mean_2_f]), power)

        m_filter = (m_1 & m_2)
        
    return m_filter

@ribo_deprecated
def get_variance_filter(kl_df, condition_1, condition_2, field, max_var=0.5):
    # first, get the field names for which we want significances
    if field == "log_translational_efficiency":
        
        # filter by both rna_abundance and ribo_abundance in both samples
        ribo_var_1_f = "ribo_abundance_var_loc_{}".format(condition_1)
        ribo_var_2_f = "ribo_abundance_var_loc_{}".format(condition_2)

        rna_var_1_f = "rna_abundance_var_loc_{}".format(condition_1)
        rna_var_2_f = "rna_abundance_var_loc_{}".format(condition_2)

        m_ribo_1 = abs(kl_df[ribo_var_1_f]) < max_var
        m_ribo_2 = abs(kl_df[ribo_var_2_f]) < max_var

        m_rna_1 = abs(kl_df[rna_var_1_f]) < max_var
        m_rna_2 = abs(kl_df[rna_var_2_f]) < max_var
        
        m_filter = (m_ribo_1 & m_ribo_2 & m_rna_1 & m_rna_2)

    else:
        var_1_f = "{}_var_loc_{}".format(field, condition_1)
        var_2_f = "{}_var_loc_{}".format(field, condition_2)

        # also get the filter
        m_1 = abs(kl_df[var_1_f]) < max_var
        m_2 = abs(kl_df[var_2_f]) < max_var

        m_filter = (m_1 & m_2)
        
    return m_filter

@ribo_deprecated
def get_mean_filter(kl_df, condition_1, condition_2, field, min_mean=1):
    # first, get the field names for which we want significances
    if field == "log_translational_efficiency":
        
        # filter by both rna_abundance and ribo_abundance in both samples
        ribo_mean_1_f = "ribo_abundance_mean_loc_{}".format(condition_1)
        ribo_mean_2_f = "ribo_abundance_mean_loc_{}".format(condition_2)

        rna_mean_1_f = "rna_abundance_mean_loc_{}".format(condition_1)
        rna_mean_2_f = "rna_abundance_mean_loc_{}".format(condition_2)

        m_ribo_1 = abs(kl_df[ribo_mean_1_f]) > min_mean
        m_ribo_2 = abs(kl_df[ribo_mean_2_f]) > min_mean

        m_rna_1 = abs(kl_df[rna_mean_1_f]) > min_mean
        m_rna_2 = abs(kl_df[rna_mean_2_f]) > min_mean
        
        m_filter = (m_ribo_1 & m_ribo_2 & m_rna_1 & m_rna_2)

    else:
        mean_1_f = "{}_mean_loc_{}".format(field, condition_1)
        mean_2_f = "{}_mean_loc_{}".format(field, condition_2)

        # also get the filter
        m_1 = abs(kl_df[mean_1_f]) > min_mean
        m_2 = abs(kl_df[mean_2_f]) > min_mean

        m_filter = (m_1 & m_2)
        
    return m_filter

@ribo_deprecated
def get_random_kl_divergence(kl_df, mean_1_f, scale_1_f, mean_2_f, scale_2_f, strategy='sampling'):
    import numpy as np
    import scipy.stats
    import misc.math_utils as math_utils

    if strategy == 'filtering':
        m_filter = [False] * len(kl_df)

        while sum(m_filter) == 0:
            x = np.random.randint(len(kl_df))
            row = kl_df.iloc[x]
        
            mean_1 = row[mean_1_f]
            scale_1 = row[scale_1_f]
            p = (mean_1, scale_1)
        
            mean_2 = row[mean_2_f]
            scale_2 = row[scale_2_f]

            m_min_scale = kl_df[scale_2_f] > 0.5*scale_2
            m_max_scale = kl_df[scale_2_f] < 2*scale_2
            m_scale = m_min_scale & m_max_scale

            m_min_mean = kl_df[mean_2_f] > 0.5*mean_2
            m_max_mean = kl_df[mean_2_f] < 2*mean_2
            m_mean = m_min_mean & m_max_mean

            m_filter = m_mean & m_scale

        indices = np.where(m_filter)[0]
        y = np.random.choice(indices)
        
        #y = np.random.randint(len(kl_df))
        row = kl_df.iloc[y]
        
        mean_2 = row[mean_2_f]
        scale_2 = row[scale_2_f]
        q = (mean_2, scale_2)

    elif strategy == 'sampling':
        x = np.random.randint(len(kl_df))
        row = kl_df.iloc[x]

        mean_1 = row[mean_1_f]
        scale_1 = row[scale_1_f]
        p = (mean_1, scale_1)
    
        mean_2 = row[mean_2_f]
        scale_2 = row[scale_2_f]

        means = kl_df[mean_2_f]

        # we take the sqrt because scipy uses std, but we use var
        #unnormalized_likelihoods = scipy.stats.norm.pdf(means, loc=mean_1, scale=np.sqrt(scale_1))
        #unnormalized_likelihoods = scipy.stats.cauchy.pdf(means, loc=mean_1, scale=np.sqrt(scale_1))

        # df=1 is the same as a cauchy
        df = 1
        unnormalized_likelihoods = scipy.stats.t.pdf(means, df, loc=mean_1, scale=np.sqrt(scale_1))
        normalized_likelihoods = unnormalized_likelihoods / np.sum(unnormalized_likelihoods)
        y = np.random.choice(len(normalized_likelihoods), p=normalized_likelihoods)
        
        row = kl_df.iloc[y]
        
        mean_2 = row[mean_2_f]
        scale_2 = row[scale_2_f]
        q = (mean_2, scale_2)

    elif strategy == "random":
        x = np.random.randint(len(kl_df))
        row = kl_df.iloc[x]

        mean_1 = row[mean_1_f]
        scale_1 = row[scale_1_f]
        p = (mean_1, scale_1)
        
        y = np.random.randint(len(kl_df))
        row = kl_df.iloc[y]
        
        mean_2 = row[mean_2_f]
        scale_2 = row[scale_2_f]
        q = (mean_2, scale_2)

    else:
        msg = "Unrecognized permutation test strategy: {}".format(strategy)
        raise ValueError(msg)





    kl = math_utils.calculate_symmetric_kl_divergence(p, q, math_utils.calculate_univariate_gaussian_kl)

    return kl

@ribo_deprecated
def get_background_kl_distribution(batch, filtered_kl_df, condition_1, condition_2, field,
                                   num_random_samples=10000, seed=8675309, use_progress_bar=False):
    
    import numpy as np
    import tqdm

    if seed is not None:
        np.random.seed(seed)

    random_kls = []
        
    # first, get the field names for which we want significances
    if field == "log_translational_efficiency":
        mean_1_f = "{}_loc_{}".format(field, condition_1)
        scale_1_f = "{}_scale_{}".format(field, condition_1)

        mean_2_f = "{}_loc_{}".format(field, condition_2)
        scale_2_f = "{}_scale_{}".format(field, condition_2)

    else:
        mean_1_f = "{}_mean_loc_{}".format(field, condition_1)
        scale_1_f = "{}_var_loc_{}".format(field, condition_1)

        mean_2_f = "{}_mean_loc_{}".format(field, condition_2)
        scale_2_f = "{}_var_loc_{}".format(field, condition_2)
    
    if use_progress_bar:
        iter_range = tqdm.trange(num_random_samples)
    else:
        iter_range = np.arange(num_random_samples)

    for i in iter_range:
        kl = get_random_kl_divergence(filtered_kl_df, mean_1_f, scale_1_f, mean_2_f, scale_2_f)
        random_kls.append(kl)
                
    return random_kls


@ribo_deprecated
def get_transcript_pvalues(kl_df, condition_1, condition_2, field, 
                min_mean=1, max_var=None, var_power=None,
                num_random_samples=10000, seed=8675309, num_cpus=1, num_groups=500):
    
    import numpy as np
    import misc.parallel as parallel
    import misc.utils as utils

    np.random.seed(seed)

    m_mean_filter = get_mean_filter(kl_df, condition_1, condition_2, 
            field, min_mean=min_mean)

    m_var_filter = True
    if max_var is not None:
        m_var_filter = get_variance_filter(kl_df, condition_1, condition_2, 
            field, max_var=max_var)
    
    m_var_power_filter = True
    if var_power is not None:
        m_var_power_filter = get_variance_power_filter(kl_df, condition_1, condition_2, 
            field, power=var_power)

    m_filter = m_mean_filter & m_var_filter & m_var_power_filter

    msg = "Total transcripts: {}. Use for sampling: {}".format(len(kl_df), sum(m_filter))
    logger.debug(msg)

    samples_per_group = np.ceil(num_random_samples / num_groups)

    # We do not need to use a seed for each group; otherwise, they all end up sampling
    # exactly the same thing.
    group_seed = None
    it = np.arange(num_cpus)
    random_kls = parallel.apply_parallel_iter(
                it,
                num_cpus,
                get_background_kl_distribution, 
                kl_df[m_filter],
                condition_1, condition_2, field, samples_per_group, group_seed,
                progress_bar=True, num_groups=num_groups)
   
    random_kls = utils.flatten_lists(random_kls)
    kls = np.array(sorted(random_kls))

    kl_field_name = "{}_{}_{}_kl_divergence".format(field, condition_1, condition_2)
    kl_field = kl_df[kl_field_name]

    pvals = kl_field.apply(get_pvalue, args=(kls,))
    
    return m_filter, pvals

@ribo_deprecated
def get_significant_differences(condition_1, condition_2, pval_df, 
                                alpha=0.05, min_rpkm_mean=None, max_rpkm_var=None,var_power=None):

    """ This function extracts the transcripts from pval_df which are
        significantly differentially "expressed" between the two given
        conditions (see below for the considered types of "expression").

        The function first filters the pval list to ensure the specified
        thresholds are met (min_rpkm_mean, max_rpkm_var, var_power). It
        then extracts the transcripts which have the specified significance
        level (alpha) or better (less) for log_transclational_efficiency,
        rna_abundance, ribo_abundance. Finally, the function returns each of
        the filters as boolean arrays.

        This function is meant to be used with the output of the 
        estimate-kl-pvalues script from the ribo-te package.

        This script uses a permutation test approach; therefore, multiple test
        correction of the pvalues *is not* required.

        Args:
            condition_1, condition_2 (strings): the name of the conditions

            pval_df (pd.DataFrame): a dataframe, which is just the output of
                the estimate-kl-pvalues script

            alpha (float): the significance value for filtering

            min_rpkm_mean, max_rpkm_var, var_power (floats): the values for filtering,
                or None if the relevant filter should not be applied.

        Returns:
            All of the return values are boolean masks of pval_df.

            m_te_filter: the transcripts which meet the filters for both RNA-seq
                and riboseq

            m_rna_filter: the transcripts which meet the filter for RNA-seq (they
                may or may not meet the riboseq filter)

            m_ribo_filter: the transcripts which meet the filter for riboseq (they
                may or may not meet the RNA-seq filter)

            m_te_sig: the transcripts which meet m_te_filter and have a significant
                KL-divergence (according to the pvalues) for log_translational_efficiency

            m_rna_sig: the transcripts which meet m_rna_filter and have a significant
                KL-divergence (according to the pvalues) for rna_abundance

            m_ribo_sig: the transcripts which meet m_ribo_filter and have a significant
                KL-divergence (according to the pvalues) for ribo_abundance

        Imports:
            numpy

    """
    import numpy as np
    
    te_kl_field = "log_translational_efficiency_{}_{}_kl_divergence".format(condition_1, condition_2)
    
    kl = pval_df[te_kl_field]

    if min_rpkm_mean is not None:
        field = "log_translational_efficiency"
        m_te_mean_filter = get_mean_filter(pval_df, condition_1, condition_2, 
            field, min_mean=min_rpkm_mean)
        
        field = "rna_abundance"
        m_rna_mean_filter = get_mean_filter(pval_df, condition_1, condition_2, 
            field, min_mean=min_rpkm_mean)
        
        field = "ribo_abundance"
        m_ribo_mean_filter = get_mean_filter(pval_df, condition_1, condition_2, 
            field, min_mean=min_rpkm_mean)
    else:
        m_te_mean_filter = True
        m_rna_mean_filter = True
        m_ribo_mean_filter = True
        
    if max_rpkm_var is not None:
        field = "log_translational_efficiency"
        m_te_var_filter = get_variance_filter(pval_df, condition_1, condition_2, 
            field, max_var=max_rpkm_var)
        
        field = "rna_abundance"
        m_rna_var_filter = get_variance_filter(pval_df, condition_1, condition_2, field, 
            max_var=max_rpkm_var)
        
        field = "ribo_abundance"
        m_ribo_var_filter = get_variance_filter(pval_df, condition_1, condition_2, field, 
            max_var=max_rpkm_var)
    else:
        m_te_var_filter = True
        m_rna_var_filter = True
        m_ribo_var_filter = True
        
    if var_power is not None:
        field = "log_translational_efficiency"
        m_te_var_power_filter = get_variance_power_filter(pval_df, condition_1, condition_2, 
            field, power=var_power)
        
        field = "rna_abundance"
        m_rna_var_power_filter = get_variance_power_filter(pval_df, condition_1, condition_2, 
            field, power=var_power)
        
        field = "ribo_abundance"
        m_ribo_var_power_filter = get_variance_power_filter(pval_df, condition_1, condition_2, 
            field, power=var_power)
    else:
        m_te_var_power_filter = True
        m_rna_var_power_filter = True
        m_ribo_var_power_filter = True
        
    field = "log_translational_efficiency"
    te_pval_field = "{}_{}_{}_pvalue".format(field, condition_1, condition_2)
    
    field = "rna_abundance"
    rna_pval_field = "{}_{}_{}_pvalue".format(field, condition_1, condition_2)
    
    field = "ribo_abundance"
    ribo_pval_field = "{}_{}_{}_pvalue".format(field, condition_1, condition_2)
    
    m_te_filter = m_te_mean_filter & m_te_var_filter & m_te_var_power_filter
    m_rna_filter = m_rna_mean_filter & m_rna_var_filter & m_rna_var_power_filter
    m_ribo_filter = m_ribo_mean_filter & m_ribo_var_filter & m_ribo_var_power_filter

    m_te_sig = (pval_df[te_pval_field] < alpha) & m_te_filter
    m_rna_sig = (pval_df[rna_pval_field] < alpha) & m_rna_filter
    m_ribo_sig = (pval_df[ribo_pval_field] < alpha) & m_ribo_filter
    
    filters = (m_te_filter, m_rna_filter, m_ribo_filter, m_te_sig, m_rna_sig, m_ribo_sig)

    filters= [ np.array(f) for f in filters ]
    return filters

@ribo_deprecated
def get_significance_filter(filters, field, significant_only=True):
    """ This function returns the appropriate mask to filter on significance
        of the given field. It assumes the filters are in the same order as the
        output of get_significant_differences.

        Parameters
        ----------
            filters : tuple
                The result of the call to get_significant_differences

            field : string
                The name of the field on which to filter. Valid options are:
                    * ribo
                    * rna
                    * te

            is_significant : bool
                Whether to return the "significant" filter (True, default) or 
                the "basic" filter

        Returns
        -------
            significant_only : boolean mask
                The appropriate mask for filtering for significance based on the
                given field.
    """

    # just map from the field to the index of the significant filters
    index_map = {
        "te": 0,
        "rna": 1,
        "ribo": 2
    }

    index = index_map[field]
    if significant_only:
        index += 3

    return filters[index]

@ribo_deprecated
def get_up_and_down_masks(condition_1, condition_2, pval_df):
    """ This function finds all of the transcripts which are, respectively
        higher or lower in the first condition. That is, "up" and "down"
        are respective to condition_1.

        This function is meant to be used with the output of the 
        estimate-kl-pvalues script from the ribo-te package.

        Args:
            condition_1, condition_2 (strings): the name of the conditions

            pval_df (pd.DataFrame): a dataframe, which is just the output of
                the estimate-kl-pvalues script
                        
        Returns:
            All of the return values are boolean masks of pval_df.

            m_te_up, m_te_down: The transcripts which have higher or lower TE
                in the first condition, respectively.

            m_rna_up, m_rna_down: The transcripts which have higher or lower
                RNA-seq RPKM in the first condition, respectively.

            m_ribo_up, m_ribo_down: The transcripts which have higher or lower
                riboseq RPKM in the first condition, respectively.

    """
    import numpy as np

    te_1 = 'log_translational_efficiency_loc_{}'.format(condition_1)
    te_2 = 'log_translational_efficiency_loc_{}'.format(condition_2)
    
    rna_1 = 'rna_abundance_mean_loc_{}'.format(condition_1)
    rna_2 = 'rna_abundance_mean_loc_{}'.format(condition_2)    
    
    ribo_1 = 'ribo_abundance_mean_loc_{}'.format(condition_1)
    ribo_2 = 'ribo_abundance_mean_loc_{}'.format(condition_2)
    
    m_te_up = pval_df[te_1] > pval_df[te_2]
    m_te_down = ~m_te_up
    
    m_rna_up = pval_df[rna_1] > pval_df[rna_2]
    m_rna_down = ~m_rna_up
    
    m_ribo_up = pval_df[ribo_1] > pval_df[ribo_2]
    m_ribo_down = ~m_ribo_up
    
    up_down_masks = m_te_up, m_te_down, m_rna_up, m_rna_down, m_ribo_up, m_ribo_down
    up_down_masks = [ np.array(f) for f in up_down_masks ]
    return up_down_masks

@ribo_deprecated
def get_up_down_filter(filters, field, direction):
    """ This function returns the appropriate mask to filter on the given field
        in the given direction. It assumes the filters are in the same order as
        the output of get_up_and_down_masks.

        Parameters
        ----------
            filters : tuple
                The result of the call to get_up_and_down_masks

            field : string
                The name of the field on which to filter. Valid options are:
                    * ribo
                    * rna
                    * te

            direction : string
                The direction in which to filter. Valid options are:
                    * up
                    * down

        Returns
        -------
            significance_mask : boolean mask
                The appropriate mask for filtering for significance based on the
                given field.
    """

    # just map from the field to the index of the significant filters
    field_map = {
        "te": 0,
        "rna": 2,
        "ribo": 4
    }

    direction_map = {
        "up": 0,
        "down": 1
    }

    index = field_map[field] + direction_map[direction]

    return filters[index]

def melt_te_df(te):
    """ Melt a data frame from the translational efficiency estimations
    to a long df suitable for use with seaborn, etc.
    """
    
    # we only want to keep the mean and var estimates
    mean_fields_to_keep = [mean_field_map[f] for f in fields]
    var_fields_to_keep = [var_field_map[f] for f in fields]
    fields_to_keep = mean_fields_to_keep + var_fields_to_keep

    # but we need this as a hierarchical index
    mean_fields = [(f, 'mean') for f in fields]
    var_fields = [(f, 'var') for f in fields]
    hierarchical_fields = mean_fields + var_fields

    # drop the rest of the fields, except gene_id
    te_df = te.set_index('gene_id')
    te_df = te_df[fields_to_keep]

    # add the multi-index for the columns
    te_df.columns = pd.MultiIndex.from_tuples(hierarchical_fields)

    # bring the gene_id back
    te_df = te_df.stack(level=0)
    te_df.index.names = ["gene_id", "field"]
    te_df = te_df.reset_index(drop=False)

    # go ahead and add the pretty name
    te_df['field_name'] = te_df['field'].map(field_map)

    return te_df

def get_bitseq_estimates(
        config,
        isoform_strategy,
        bitseq_id_field='transcript_id',
        strings_to_remove=['.cds-only', '.merged']):
    """ Load the bitseq abundance estimates into a single long data frame.

    Parameters
    ----------
    config: dict-like
        The configuration for the project, presumably from the yaml file

    isoform_strategy: str
        The strategy for handling transcript isoforms

    bitseq_id_field: str
        Name for the "transcript_id" field (second column) in bitseq tr file

    strings_to_remove: list of strings
        A list of strings to replace with "" in the bitseq ids

    Returns
    -------
    bitseq_estimates: pd.DataFrame
        A data frame containing the following columns

            * rpkm_{mean,var}: the bitseq estimates
            * sample: the name of the respective sample
            * type: "ribo" or "rna"
    """
    import misc.bio as bio
    import tqdm
    
    msg = "Reading the bitseq tr info file"
    logger.info(msg)

    # check which transcript file to load
    is_merged = False
    if isoform_strategy == "merged":
        is_merged = True

    # and get the file
    transcript_fasta = filenames.get_transcript_fasta(
        config['genome_base_path'], 
        config['genome_name'], 
        is_annotated=True, 
        is_merged=is_merged, 
        is_cds_only=True
    )

    tr_info = filenames.get_bitseq_transcript_info(transcript_fasta)
    bitseq_tr = bio.read_bitseq_tr_file(tr_info)

    # we need to remove all of the indicated strings from the ids
    for to_remove in strings_to_remove:
        tids = bitseq_tr['transcript_id'].str.replace(to_remove, "")
        bitseq_tr['transcript_id'] = tids

    bitseq_tr = bitseq_tr.rename(columns={'transcript_id': bitseq_id_field})

    note = config.get('note', None)
    
    all_dfs = []
    
    msg = "Reading riboseq BitSeq estimates"
    logger.info(msg)
    
    is_unique = 'keep_riboseq_multimappers' not in config
    
    it = tqdm.tqdm(config['riboseq_samples'].items())
    for name, file in it:

        lengths, offsets = get_periodic_lengths_and_offsets(
            config, 
            name,
            isoform_strategy=isoform_strategy, 
            is_unique=is_unique
        )

        bitseq_rpkm_mean = filenames.get_riboseq_bitseq_rpkm_mean(
            config['riboseq_data'], 
            name, 
            is_unique=is_unique, 
            is_transcriptome=True, 
            is_cds_only=True,
            length=lengths, 
            offset=offsets, 
            isoform_strategy=isoform_strategy, 
            note=note
        )

        field_names = ['rpkm_mean', 'rpkm_var']
        bitseq_rpkm_mean_df = bio.read_bitseq_means(
            bitseq_rpkm_mean,
            names=field_names
        )

        bitseq_rpkm_mean_df['sample'] = name
        bitseq_rpkm_mean_df['type'] = 'ribo'
        bitseq_rpkm_mean_df[bitseq_id_field] = bitseq_tr[bitseq_id_field]

        all_dfs.append(bitseq_rpkm_mean_df)


    # now, the rnaseq    
    msg = "Reading RNA-seq BitSeq estimates"
    logger.info(msg)
    
    is_unique = ('remove_rnaseq_multimappers' in config)
    
    it = tqdm.tqdm(config['rnaseq_samples'].items())
    for name, data in it:

        bitseq_rpkm_mean = filenames.get_rnaseq_bitseq_rpkm_mean(
            config['rnaseq_data'], 
            name, 
            is_unique=is_unique, 
            is_transcriptome=True, 
            is_cds_only=True,
            isoform_strategy=isoform_strategy, 
            note=note
        )

        field_names = ['rpkm_mean', 'rpkm_var']
        bitseq_rpkm_mean_df = bio.read_bitseq_means(
            bitseq_rpkm_mean,
            names=field_names
        )

        bitseq_rpkm_mean_df['sample'] = name
        bitseq_rpkm_mean_df['type'] = 'rna'
        bitseq_rpkm_mean_df[bitseq_id_field] = bitseq_tr[bitseq_id_field]

        all_dfs.append(bitseq_rpkm_mean_df)

    
    msg = "Joining estimates into long data frame"
    logger.info(msg)
    
    long_df = pd.concat(all_dfs)
    long_df = long_df.reset_index(drop=True)
    return long_df

def update_gene_id_from_transcript_id(df:pd.DataFrame, config:dict, args=None):
    """ Assuming "gene_id" is actually a transcript id, replace it
    with the actual gene identifier.
    
    This function is used in the case of the "all" isoform
    strategy when downstream analysis actually needs a gene
    identifier.
    
    Parameters
    ----------
    df: pd.DataFrame
        A data frame which contains a "gene_id" field which actually
        contains transcript identifiers. For example, the latter parts
        of the B-tea pipeline produce data frames like this with
        the "all" isoform strategy
        
    config: dict
        Configuration options
        
    args: argparse.Namespace or None
        The logging options from the command line. pyensembl likes
        to overwrite these, so they will be reset.
        
    Returns
    -------
    updated_df: pd.DataFrame
        A data frame in which the 'gene_id' column is moved to a
        'transcript_id' column, and the 'gene_id' column is updated
        to include  actual gene identifiers
    """
    import misc.bio_utils.pyensembl_utils as pyensembl_utils

    msg = "Loading Ensembl annotations"
    logger.info(msg)
    
    ensembl = pyensembl_utils.get_genome(
        config['genome_name'],
        config['gtf'],
        logging_args = args
    )
    
    msg = "Finding the gene ids for each transcript id"
    logger.info(msg)
    
    gene_ids = set(df['gene_id'])
    transcript_gene_mapping = pyensembl_utils.get_gene_ids_of_transcript_ids(
        gene_ids, ensembl)
    
    msg = "Adding gene ids to data frame"
    logger.info(msg)
    
    df['transcript_id'] = df['gene_id']
    df = df.drop('gene_id', 1)
    df = df.merge(transcript_gene_mapping, on='transcript_id')
    
    return df

###
#   These are functions for retrieving the dominant isoform for
#   each gene and condition.
###

def _get_matching_condition(row, condition_field, config):
    condition = row[condition_field]
    field = row['field']
    
    # use the ribo conditions for te
    if field == "te":
        field = "ribo"
    
    return ribo_utils.get_criterion_condition(condition, field, config)

def _add_matching_conditions(pvalues, config):
    """ Add the "matching" conditions for both conditions. """

    # turn off logging; we already know we have matching conditions
    logger_level = logger.getEffectiveLevel()

    logger.setLevel("WARNING")
    matching_condition_1 = parallel.apply_df_simple(
        pvalues,
        _get_matching_condition,
        "condition_1",
        config
    )

    matching_condition_2 = parallel.apply_df_simple(
        pvalues,
        _get_matching_condition,
        "condition_2",
        config
    )

    pvalues['matching_condition_1'] = matching_condition_1
    pvalues['matching_condition_2'] = matching_condition_2

    logger.setLevel(logger_level)
    
    return pvalues

def _add_transcript_id(pvalues, abundances):
    """ Use the gene ID an dominant isoform information
    to pull back the transcript id for each "matching" condition.
    """
    
    left_on=['matching_condition_1', 'gene_id', 'field']
    right_on=['condition', 'gene_id', 'field']

    pvalues = pvalues.merge(abundances, left_on=left_on, right_on=right_on)
    pvalues = pvalues.rename(columns={"transcript_id": "transcript_id_1"})
    pvalues = pvalues.drop('condition', 1)

    left_on=['matching_condition_2', 'gene_id', 'field']
    pvalues = pvalues.merge(abundances, left_on=left_on, right_on=right_on)
    pvalues = pvalues.rename(columns={"transcript_id": "transcript_id_2"})
    pvalues = pvalues.drop('condition', 1)
    
    return pvalues

def get_dominant_transcript_ids(pvalues:pd.DataFrame, config:dict):
    """ Add the transcript id for the dominant isoform in each condition.
    
    This function is really only intended to be used with the final pvalues
    data frame from B-tea.
    """
    
    # now, we need to get the transcript ids for condition_1 and condition_2
    abundance_fields_to_keep = [
        'type',
        'transcript_id',
        'gene_id',
        'condition'
    ]

    msg = "Reading abundances"
    logger.info(msg)

    abundances = filenames.get_abundances(
        config['translational_efficiency_data'],
        isoform_strategy=args.isoform_strategy,
        note=note
    )

    abundances = pd.read_csv(abundances)
    abundances = abundances[abundance_fields_to_keep]
    abundances = abundances.drop_duplicates()
    abundances = abundances.rename(columns={"type": "field"})
    
    pvalues = _add_matching_conditions(pvalues, config)
    pvalues = _add_transcript_id(pvalues, abundances)
    return pvalues

###
#   End of dominant isoform extraction functions
###


def get_overlap_data_frame(unique_file, multimappers_file):
    import pandas as pd
    import misc.bio_utils.bed_utils as bed_utils

    msg = "Reading predictions with unique mappers"
    logger.info(msg)
    unique = bed_utils.read_bed(unique_file)

    msg = "Reading predictions with multimappers"
    logger.info(msg)
    multimappers = bed_utils.read_bed(multimappers_file)

    msg = "Splitting predictions with multimappers"
    logger.info(msg)
    multimappers_exons = bed_utils.split_bed12(multimappers)

    msg = "Splitting predictions with unique mappers"
    logger.info(msg)
    unique_exons = bed_utils.split_bed12(unique)

    msg = "Finding overlap"
    logger.info(msg)
    overlap = bed_utils.get_bed_overlaps(unique, multimappers, 
        exons_a=unique_exons, exons_b=multimappers_exons)

    msg = "Constructing data frame with overlaps and ORFs from each prediction set"
    logger.info(msg)

    unique_with_overlap = {o.a_info for o in overlap}
    multimapper_with_overlap = {o.b_info for o in overlap}

    overlap_df = pd.DataFrame(overlap)
    overlap_df = overlap_df.rename(columns={
        "a_fraction":"Unique Coverage", "b_fraction": "Multimapper Coverage"
    })
    overlap_df['category'] = 'overlap'

    m_unique_with_overlap = unique['id'].isin(unique_with_overlap)
    m_multimapper_with_overlap = multimappers['id'].isin(multimapper_with_overlap)

    unique_no_overlap = unique.loc[~m_unique_with_overlap, 'id']
    multimapper_no_overlap = multimappers.loc[~m_multimapper_with_overlap, 'id']

    unique_df = pd.DataFrame()
    unique_df['a_info'] = unique_no_overlap
    unique_df['b_info'] = ""
    unique_df['overlap'] = 0
    unique_df['Unique Coverage'] = 0
    unique_df['Multimapper Coverage'] = 0
    unique_df['category'] = 'unique_only'

    multimapper_df = pd.DataFrame()
    multimapper_df['a_info'] = ""
    multimapper_df['b_info'] = multimapper_no_overlap
    multimapper_df['overlap'] = 0
    multimapper_df['Unique Coverage'] = 0
    multimapper_df['Multimapper Coverage'] = 0
    multimapper_df['category'] = 'multimapper_only'
    joined_df = pd.concat([overlap_df, unique_df, multimapper_df])

    msg = "Adding expression, etc., to data frame"
    logger.info(msg)

    joined_df = joined_df.merge(unique, left_on="a_info", right_on="id", 
        suffixes=['', '_unique'], how='left')

    to_rename = {c: "{}_unique".format(c) for c in unique.columns}
    joined_df = joined_df.rename(columns=to_rename)

    joined_df = joined_df.merge(multimappers, left_on="b_info", right_on="id", 
        suffixes=['', '_multimappers'], how='left')

    to_rename = {c: "{}_multimappers".format(c) for c in multimappers.columns}
    joined_df = joined_df.rename(columns=to_rename)

    ui = joined_df['x_1_sum_unique'].divide(joined_df['profile_sum_unique']) 
    joined_df['inframe_unique'] = ui

    mi = joined_df['x_1_sum_multimappers'].divide(joined_df['profile_sum_multimappers'])
    joined_df['inframe_multimappers'] = mi

    joined_df = joined_df.fillna(0)

    return joined_df
