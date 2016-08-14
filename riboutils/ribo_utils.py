import logging
import os
import pandas as pd

import riboutils.ribo_filenames as filenames

logger = logging.getLogger(__name__)

###
#   The following labels are used to group similar ORF types.
###
orf_types = ['canonical', 'canonical_extended', 'canonical_truncated', 
    'five_prime', 'three_prime', 'noncoding', 'novel', 'five_prime_overlap',
    'suspect_overlap', 'three_prime_overlap', 'within']


orf_type_labels = ['canonical', 'canonical_variant', 'five_prime', 
    'three_prime', 'noncoding', 'other', 'novel']

orf_type_labels_mapping = {
    'canonical': ['canonical'],
    'canonical_variant': ['canonical_extended', 'canonical_truncated'],
    'five_prime': ['five_prime'],
    'three_prime': ['three_prime'],
    'noncoding': ['noncoding'],
    'novel': ['novel'],
    'other': ['five_prime_overlap', 'suspect_overlap', 'three_prime_overlap', 'within']
}

###
#   The following functions are all used for parsing replicates, etc., from the config file.
###

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

def get_periodic_lengths_and_offsets(config, name, do_not_call=False, is_merged=False):
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

        Args:
            config (dictionary) : the configuration information(see description)

            name (string) : the name of the dataset in question

            do_not_call (bool) : whether the metagene bf file should exist. If false,
                then dummy values are returned (and a warning message is printed).

            is_merged (bool) : whether the "merged" transcripts are used (i.e., is
                this for rpbp or ribo-te)

        Returns:
            lengths (list of strings) : all of the periodic read lengths

            offsets (list of strings) : the corresponding P-site offsets for the
                read lengths. 

        Imports:
            numpy
            scipy.stats
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

    periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], name, 
        is_unique=True, is_merged=is_merged, note=note_str)
    
    if not os.path.exists(periodic_offsets):
        msg = ("The periodic offsets file does not exist. Please ensure the select-periodic-offsets "
            "script completed successfully or specify the \"use_fixed_lengths\", \"lengths\", and "
            "\"offsets\" values in the configuration file. '{}'".format(periodic_offsets))

        if do_not_call:
            msg = msg +  ("\nThe --do-not-call flag was given, so \"dummy\" default lengths will be "
                "used to check the remaining calls.\n")

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

        msg = "Num nans: {}, num predictions: {}".format(num_nans, num_predictions)
        logger.debug(msg)

        msg = ("Using the likelihood filter. min_mean: {}, min_likelihood: {}"
            .format(min_bf_mean, min_bf_likelihood))
        logger.debug(msg)

        max_likelihood = max(likelihood[~nans])
        msg = "Maximum likelihood: {}".format(max_likelihood)
        logger.debug(msg)

        # now filter
        m_bf_likelihood = likelihood > min_bf_likelihood

    if (max_bf_var is None) and (min_bf_likelihood is None):
        m_bf_mean = bf['highest_peak_bf_mean'] > min_bf_mean

    filtered_periodic_offsets = offsets_df[m_count & m_bf_mean & m_bf_var & m_bf_likelihood]

    offsets = filtered_periodic_offsets['highest_peak_offset']
    lengths = filtered_periodic_offsets['length']

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

    for i, a in enumerate(tqdm.tqdm(alignments, leave=True, file=sys.stdout, total=num_alignments)):
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

        likelihood = 1-scipy.stats.norm.cdf(min_bf_mean, bf['bayes_factor_mean'], np.sqrt(bf['bayes_factor_var']))

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
                            chisq_alpha=default_chisq_alpha):
    """ This function applies a set of filters to ORFs to select those which
        are predicted as "translated." This function selects translated ORFs
        based on the Bayes factor estimates and the chi-square p-values. ORFs
        must pass all of the relevant features to be selected as "translated."
        Finally, among all ORFs which share a stop codon, only longest
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

        Returns:
            longest_orfs (pd.DataFrame) : all longest ORFs which meet the profile,
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
    import numpy as np
    import scipy.stats

    msg = "Finding all longest ORFs with signal"
    logger.info(msg)

    m_base = get_base_filter(bf, min_signal, min_length)

    longest_orfs = bio.get_longest_features_by_end(bf[m_base])
    
    # create the selected ORFs
    m_bf =  get_bf_filter(bf, min_bf_mean, max_bf_var, min_bf_likelihood)

    # apply all the filters
    m_bf_predicted = m_base & m_bf

    bf_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_bf_predicted])

    M = len(longest_orfs)
    # for the bonferroni correction, we only correct for the number of tests we actually consider
    # that is, we only correct for orfs which pass the base filter
    corrected_significance_level = chisq_alpha / M

    msg = "Corrected significance level: {}".format(corrected_significance_level)
    logger.debug(msg)
    
    m_chisq_pval = bf['chi_square_p'] < corrected_significance_level
    m_chisq_predicted = m_base & m_chisq_pval

    chisq_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_chisq_predicted])
    
    return (longest_orfs, bf_longest_predicted_orfs, chisq_longest_predicted_orfs)
    

###
# The following functions are all related. They are used to estimate p-values
# for the KL-divergence values calculated for translational efficiency (only).
###
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

    return kl, p, q

def get_background_kl_distribution(batch, filtered_kl_df, condition_1, condition_2, field,
                                   num_random_samples=10000, seed=8675309, use_progress_bar=False):
    
    import numpy as np
    import tqdm

    if seed is not None:
        np.random.seed(seed)

    random_kls = []
    random_ps = []
    random_qs = []
        
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
        kl, p, q = get_random_kl_divergence(filtered_kl_df, mean_1_f, scale_1_f, mean_2_f, scale_2_f)
        random_kls.append(kl)
        random_ps.append(p)
        random_qs.append(q)
        
    return random_kls, random_ps, random_qs

def get_pvalue(val, kls):
    import numpy as np

    pos = np.searchsorted(kls, val)
    p = 1 - (pos/len(kls))
    return p

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
   
    print(len(random_kls))
    random_ks = utils.flatten_lists(random_kls)
    print(len(random_ks))
    #random_ks = np.array(random_kls, dtype=object)
    #print(random_ks.shape)
    random_kls = random_ks[0::3]
    random_ps = random_ks[1::3]
    random_qs = random_ks[2::3]

    random_kls = np.concatenate(random_kls)
    random_ps = np.concatenate(random_ps)
    random_qs = np.concatenate(random_qs)

    print(random_kls.shape)
    print(random_qs.shape)
    print(random_ps.shape)
    #print(random_kls)
    #print(random_ps)
    #print(random_qs)
    kls = np.array(sorted(random_kls))

    kl_field_name = "{}_{}_{}_kl_divergence".format(field, condition_1, condition_2)
    kl_field = kl_df[kl_field_name]

    pvals = kl_field.apply(get_pvalue, args=(kls,))
    
    return m_filter, pvals, random_kls, random_ps.tolist(), random_qs.tolist()

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
    
    return (m_te_filter, m_rna_filter, m_ribo_filter, m_te_sig, m_rna_sig, m_ribo_sig)

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
    
    return m_te_up, m_te_down, m_rna_up, m_rna_down, m_ribo_up, m_ribo_down

###
#   These functions are used for estimating significance (p-values) of the overlap
#   with Mackowiak ORFs via a permutation test.
#   
###

def get_random_mackowiak_exact_matches(bf_filtered, bf_exons, num_random_choices, sorfs_exons):
    """ This function performs one step of a permutation test for determining
        p-values for the overlap of Mackowiak ORFs. In particular, given a list
        of ORFs and a number of random ORFs to choose, this function selects 
        that many ORFs. It then joins them to find the "longest"
        predicted ORFs. Finally, it intersects that list with the Mackowiak ORFs.
        The function returns the number of exact matches between the randomly
        drawn ORFs and the Mackowiak ORFs.
        
        Args:
            bf_filtered (pd.DataFrame): a data frame of ORFs (i.e., the output of 
                estimate-orf-bayes-factors). This should already include
                basic filtering (length, etc.)
                
            num_random_choices (int): the number of random ORFs to choose (this
                should be the same as the number of ORFs predicted as translated
                before merging into the "longest" ORFs)
                
            sorfs_exons (pd.DataFrame): the split exons of the Mackowiak ORFs
                        
        Returns:
            int: the number of exact matches between the randomly chosen ORFs
                and the Mackowiak ORFs
                
        Imports:
            numpy
            misc.bio_utils.bed_utils
            time
    """
    import time
    import numpy as np
    import misc.bio_utils.bed_utils as bed_utils

    start = time.perf_counter()
    
    # randomly choose some of the ORFs which passed the base filter
    random_choices = np.random.choice(len(bf_filtered), num_random_choices, replace=False)
    
    # convert the indices back into a mask we can use
    m_random = np.zeros(len(bf_filtered), dtype=bool)
    m_random[random_choices] = True

    t = time.perf_counter() - start
    msg = "Time for sampling: {}".format(t)
    logger.debug(msg)
    
    # merge the selected ORFs
    random_orfs = bed_utils.get_longest_features_by_end(bf_filtered[m_random])

    t = time.perf_counter() - start
    msg = "Time for finding longest ORFs: {}".format(t)
    logger.debug(msg)

    #random_orfs_exons = bed_utils.split_bed12(random_orfs)
    random_orfs_fields = ['id']
    random_orfs_exons = random_orfs[random_orfs_fields].merge(bf_exons, on='id')

    t = time.perf_counter() - start
    msg = "Time for splitting ORFs: {}".format(t)
    logger.debug(msg)

    # now, take the intersection with the Mackowiak ORFs
    i_random_sorfs = bed_utils.get_all_exact_bed_matches(random_orfs_exons, sorfs_exons)

    t = time.perf_counter() - start
    msg = "Time for finding exact matches: {}".format(t)
    logger.debug(msg)
    
    # and return the number of exact matches
    num_exact_matches = len(i_random_sorfs)
    return num_exact_matches


def get_mackowiak_background(batch, num_random_samples, 
                             bf_filtered,
                             num_random_choices, 
                             sorfs_file,
                             progress_bar=False):

    """ This function draws the specified number of random samples for a 
        permutation test using get_random_mackowiak_exact_matches. It is mostly
        just a wrapper around that function that is convenient for parallel
        calls.

        Args:
            batch (int): the index of this set of samples

            num_random_samples (int): the number of random samples to draw

            progress_bar (bool): whether to show a progress bar

            the rest are the same as for get_random_mackowiak_exact_matches

        Returns:
            np.array of ints: the number of exact matches from each random
                sample

        Imports:
            misc.bio_utils.bed_utils
            numpy
            tqdm (if a progress bar is used)
            sys (if a progress bar is used)
    """
    import sys
    import time

    import numpy as np

    import tqdm
    import misc.bio_utils.bed_utils as bed_utils

    start = time.perf_counter()
    
    sorfs = bed_utils.read_bed(sorfs_file)
    t = time.perf_counter() - start
    msg = "Time to read Mackowiak ORFs: {}".format(t)
    logger.debug(msg)

    sorfs_exons = bed_utils.split_bed12(sorfs)
    t = time.perf_counter() - start
    msg = "Time to split Mackowiak ORFs: {}".format(t)
    logger.debug(msg)
    
    bf_exons = bed_utils.split_bed12(bf_filtered)
    t = time.perf_counter() - start
    msg = "Time to split ORFs: {}".format(t)
    logger.debug(msg)
    
    ret = np.zeros(num_random_samples, dtype=int)
        
    it = range(num_random_samples)
    if progress_bar:
        import sys
        import tqdm
        it = tqdm.tqdm(it, leave=True, file=sys.stdout, total=num_random_samples)
    
    for i in it:
        num_exact_matches = get_random_mackowiak_exact_matches(bf_filtered,
                                                               bf_exons,
                                                               num_random_choices, 
                                                               sorfs_exons)
        ret[i] = num_exact_matches
        
    return ret
