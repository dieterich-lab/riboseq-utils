import logging
import os

import ribo_utils.filenames as filenames

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
default_min_metagene_profile_bayes_factor_mean = 5
default_max_metagene_profile_bayes_factor_var = 5

def get_periodic_lengths_and_offsets(config, name, do_not_call=False):
    import pandas as pd

    # check if we specified to just use a fixed offset and length
    if 'use_fixed_lengths' in config:
        lengths = config['lengths']
        offsets = config['offsets']

        return (lengths, offsets)

    # filter out the lengths which do not satisfy the quality thresholds
    min_metagene_profile_count = config.get(
        "min_metagene_profile_count", default_min_metagene_profile_count)

    min_metagene_profile_bayes_factor_mean = config.get(
        "min_metagene_profile_bayes_factor_mean", default_min_metagene_profile_bayes_factor_mean)

    max_metagene_profile_bayes_factor_var = config.get(
        "max_metagene_profile_bayes_factor_var", default_max_metagene_profile_bayes_factor_var)

    note_str = config.get('note', None)

    periodic_offsets = filenames.get_periodic_offsets(config['riboseq_data'], name, 
        is_unique=True, note=note_str)
    
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
    m_count = offsets_df['highest_peak_profile_sum'] > min_metagene_profile_count
    m_bf_mean = offsets_df['highest_peak_bf_mean'] > min_metagene_profile_bayes_factor_mean
    m_bf_var = offsets_df['highest_peak_bf_var'] < max_metagene_profile_bayes_factor_var

    filtered_periodic_offsets = offsets_df[m_count & m_bf_mean & m_bf_var]

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

default_min_signal = 10
default_min_bf_mean = 5
default_max_bf_var = 2
default_min_length = 200
default_chisq_alpha = 0.01

def get_predicted_orfs(bf, min_signal=default_min_signal, min_bf_mean=default_min_bf_mean, 
                       max_bf_var=default_max_bf_var, min_length=default_min_length,
                       chisq_alpha=default_chisq_alpha):

    import misc.bio as bio

    msg = "Finding all longest ORFs with signal"
    logging.info(msg)

    mask_profile = bf['profile_sum'] > min_signal
    longest_orfs = bio.get_longest_features_by_end(bf[mask_profile])
    
    # create the selected ORFs
    m_profile_sum = bf['profile_sum'] > min_signal
    m_bf_mean = bf['bayes_factor_mean'] > min_bf_mean
    m_bf_var = bf['bayes_factor_var'] < max_bf_var
    m_x1_gt_x2 = bf['x_1_sum'] > bf['x_2_sum']
    m_x1_gt_x3 = bf['x_1_sum'] > bf['x_3_sum']
    
    m_length = bf['orf_len'] > min_length

    m_bf_predicted = m_x1_gt_x2 & m_x1_gt_x3 & m_profile_sum & m_bf_mean & m_bf_var & m_length

    bf_longest_predicted_orfs = bio.get_longest_features_by_end(bf[m_bf_predicted])

    M = len(longest_orfs)
    # for the bonferroni correction, we only correct for the number of tests we actually consider
    # that is, we only correct for orfs which pass the minimum profile filter
    corrected_significance_level = chisq_alpha / sum(m_profile_sum)

    msg = "Corrected significance level: {}".format(corrected_significance_level)
    logging.debug(msg)
    
    m_chisq_pval = bf['chi_square_p'] < corrected_significance_level
    m_chisq_predicted = m_x1_gt_x2 & m_x1_gt_x3 & m_profile_sum & m_chisq_pval & m_length

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
