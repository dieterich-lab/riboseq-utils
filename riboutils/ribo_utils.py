import logging
import os
import pandas as pd

import riboutils.ribo_filenames as filenames

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
            logging.info(msg)
            return config['riboseq_biological_replicates']
        
    msg = ("Did not find 'riboseq_biological_replicates' key in config file. "
            "Using each 'riboseq_sample' as a single-condition replicate.")
    logging.info(msg)

    # create a dictionary mapping from the sample name to asingle-element list
    ret = {
        name: [name] for name, sample in config['riboseq_samples'].items()
    }
    return ret

def get_rnaseq_replicates(config):
    if 'rnaseq_biological_replicates' in config:
        if config['rnaseq_biological_replicates'] is not None:
            msg = "Found 'rnaseq_biological_replicates' key in config file"
            logging.info(msg)

            return config['rnaseq_biological_replicates']
        
    msg = ("Did not find 'rnaseq_biological_replicates' key in config file. "
            "Using each 'rnaseq_sample' as a single-condition replicate.")
    logging.info(msg)
    
    # create a dictionary mapping from the sample name to asingle-element list
    ret = {
        name: [name] for name, sample in config['rnaseq_samples'].items()
    }
    return ret

def get_matching_conditions(config):
    if 'matching_conditions' in config:
        if config['matching_conditions'] is not None:
            msg = "Found 'matching_conditions' key in config file"
            logging.info(msg)

            return config['matching_conditions']
        
    msg = ("Did not find 'matching_conditions' key in config file. Using "
            "riboseq and rnaseq conditions (biological_replicate entries) "
            "as matching conditions.")
    logging.info(msg)
    
    # otherwise, get the replicates and match key names
    riboseq_replicates = get_riboseq_replicates(config)
    rnaseq_replicates = get_rnaseq_replicates(config)
    
    matching_conditions = {
        x: [x, x] for x in riboseq_replicates if x in rnaseq_replicates
    }
    
    return matching_conditions

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
        lengths = config['lengths']
        offsets = config['offsets']

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

            logging.warning(msg)

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
        logging.debug(msg)

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
        logging.debug(msg)

        msg = ("Using the likelihood filter. min_mean: {}, min_likelihood: {}"
            .format(min_bf_mean, min_bf_likelihood))
        logging.debug(msg)

        max_likelihood = max(likelihood[~nans])
        msg = "Maximum likelihood: {}".format(max_likelihood)
        logging.debug(msg)

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
#
# This function extracts all ORFs which count as "translated", according to the
# values in the config file.
#
###

default_min_signal = None
default_min_bf_mean = 5
default_max_bf_var = None
default_min_bf_likelihood = None
default_min_length = 0
default_chisq_alpha = 0.01

def get_predicted_orfs(bf, min_signal=default_min_signal, 
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
    logging.info(msg)

    if min_signal is None:
        m_profile = bf['profile_sum'] > 0
    else:
        m_profile = bf['profile_sum'] > min_signal

    m_length = bf['orf_len'] > min_length
    m_x1_gt_x2 = bf['x_1_sum'] > bf['x_2_sum']
    m_x1_gt_x3 = bf['x_1_sum'] > bf['x_3_sum']

    m_base = m_profile & m_length & m_x1_gt_x2 & m_x1_gt_x3

    longest_orfs = bio.get_longest_features_by_end(bf[m_base])
    
    # create the selected ORFs

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
        logging.debug(msg)

        if num_nans != num_predictions:
            max_likelihood = max(likelihood[~nans])
            msg = "Maximum likelihood: {}".format(max_likelihood)
            logging.debug(msg)

        # now filter
        m_bf_likelihood = likelihood > min_bf_likelihood

    if (max_bf_var is None) and (min_bf_likelihood is None):
        m_bf_mean = bf['bayes_factor_mean'] > min_bf_mean

    # apply all the filters
    m_bf_predicted = m_base & m_bf_mean & m_bf_var & m_bf_likelihood

    bf_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_bf_predicted])

    M = len(longest_orfs)
    # for the bonferroni correction, we only correct for the number of tests we actually consider
    # that is, we only correct for orfs which pass the base filter
    corrected_significance_level = chisq_alpha / M

    msg = "Corrected significance level: {}".format(corrected_significance_level)
    logging.debug(msg)
    
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

        m_ribo_1 = kl_df[ribo_var_1_f] < np.power(kl_df[ribo_mean_1_f], power)
        m_ribo_2 = kl_df[ribo_var_2_f] < np.power(kl_df[ribo_mean_2_f], power)
        
        m_rna_1 = kl_df[rna_var_1_f] < np.power(kl_df[rna_mean_1_f], power)
        m_rna_2 = kl_df[rna_var_2_f] < np.power(kl_df[rna_mean_2_f], power)
        
        m_filter = (m_ribo_1 & m_ribo_2 & m_rna_1 & m_rna_2)

    else:
        var_1_f = "{}_var_loc_{}".format(field, condition_1)
        var_2_f = "{}_var_loc_{}".format(field, condition_2)

        mean_1_f = "{}_mean_loc_{}".format(field, condition_1)
        mean_2_f = "{}_mean_loc_{}".format(field, condition_2)


        # also get the filter
        m_1 = kl_df[var_1_f] < np.power(kl_df[mean_1_f], power)
        m_2 = kl_df[var_2_f] < np.power(kl_df[mean_2_f], power)

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

        m_ribo_1 = kl_df[ribo_var_1_f] < max_var
        m_ribo_2 = kl_df[ribo_var_2_f] < max_var

        m_rna_1 = kl_df[rna_var_1_f] < max_var
        m_rna_2 = kl_df[rna_var_2_f] < max_var
        
        m_filter = (m_ribo_1 & m_ribo_2 & m_rna_1 & m_rna_2)

    else:
        var_1_f = "{}_var_loc_{}".format(field, condition_1)
        var_2_f = "{}_var_loc_{}".format(field, condition_2)

        # also get the filter
        m_1 = kl_df[var_1_f] < max_var
        m_2 = kl_df[var_2_f] < max_var

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

        m_ribo_1 = kl_df[ribo_mean_1_f] > min_mean
        m_ribo_2 = kl_df[ribo_mean_2_f] > min_mean

        m_rna_1 = kl_df[rna_mean_1_f] > min_mean
        m_rna_2 = kl_df[rna_mean_2_f] > min_mean
        
        m_filter = (m_ribo_1 & m_ribo_2 & m_rna_1 & m_rna_2)

    else:
        mean_1_f = "{}_mean_loc_{}".format(field, condition_1)
        mean_2_f = "{}_mean_loc_{}".format(field, condition_2)

        # also get the filter
        m_1 = kl_df[mean_1_f] > min_mean
        m_2 = kl_df[mean_2_f] > min_mean

        m_filter = (m_1 & m_2)
        
    return m_filter

def get_random_kl_divergence(kl_df, mean_1_f, scale_1_f, mean_2_f, scale_2_f):
    import numpy as np
    import misc.math_utils as math_utils

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

    kl = math_utils.calculate_symmetric_kl_divergence(p, q, math_utils.calculate_univariate_gaussian_kl)

    return kl

def get_background_kl_distribution(filtered_kl_df, condition_1, condition_2, field,
                                   num_random_samples=10000, seed=8675309, use_progress_bar=False):
    
    import numpy as np
    import tqdm

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

def get_pvalue(val, kls):
    import numpy as np

    pos = np.searchsorted(kls, val)
    p = 1 - (pos/len(kls))
    return p

def get_transcript_pvalues(kl_df, condition_1, condition_2, field, 
                min_mean=1, num_random_samples=10000, seed=8675309, num_procs=1, num_groups=500):
    
    import numpy as np
    import misc.parallel as parallel
    import misc.utils as utils

    m_filter = get_mean_filter(kl_df, condition_1, condition_2, field, min_mean=min_mean)

    samples_per_group = np.ceil(num_random_samples / num_groups)

    random_kls = parallel.apply_parallel_split(
                kl_df[m_filter],
                num_procs,
                get_background_kl_distribution, 
                condition_1, condition_2, field, samples_per_group, seed,
                progress_bar=True, num_groups=num_groups)
    
    random_kls = utils.flatten_lists(random_kls)
    kls = np.array(sorted(random_kls))

    kl_field_name = "{}_{}_{}_kl_divergence".format(field, condition_1, condition_2)
    kl_field = kl_df[kl_field_name]

    pvals = kl_field.apply(get_pvalue, args=(kls,))
    
    return m_filter, pvals

def get_significant_differences(condition_1, condition_2, pval_df, 
                                alpha=0.05, min_rpkm_mean=None, max_rpkm_var=None,var_power=None):
    
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
