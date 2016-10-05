import glob
import os

### parameterized names

def get_cds_only_string(is_cds_only):
    cds_only = ""
    if is_cds_only:
        cds_only = ".cds-only"
    return cds_only

def get_chisq_string(is_chisq):
    chisq = ""
    if is_chisq:
        chisq = ".chisq"
    return chisq

def get_fastqc_name(filename):
    """ Given the sequence or alignment filename, this function extracts the name
        used by fastqc. In particular, it removes the ".fasta[.gz]" or ".fastq[.gz]"
        or ".sam" or ".bam" ending from the filename.

        Args:
            filename: the name of the sequence or alignment file (NOT including the path)

        Returns:
            string : the fastqc name

        Imports:
            misc.utils
    """

    import misc.utils as utils

    # first, get the filename
    filename = utils.get_basename(filename)
    filename = filename.replace(".fasta", "")
    filename = filename.replace(".fastq", "")
    filename = filename.replace(".sam", "")
    filename = filename.replace(".bam", "")
    filename = filename.replace(".gz", "")

    return filename

def get_filtered_string(is_filtered):
    filtered_str = ""
    if is_filtered:
        filtered_str = ".filtered"
    return filtered_str

def get_fraction_string(fraction=None):
    fraction_str = ""
    if (fraction is not None) and  (len(str(fraction)) > 0):
        fraction_str = ".frac-{}".format(fraction)
    return fraction_str

def get_grouped_string(is_grouped):
    g = ""
    if is_grouped:
        g = ".grouped"
    return g

def get_isoforms_string(is_isoform):
    i = ""
    if is_isoform:
        i = ".isoforms"
    return i


def get_length_string(length=None):
    l = ""
    if length is not None:
        if isinstance(length, (list, tuple)):
            l = "-".join(length)
            l = ".length-{}".format(l)
        else:
            l = ".length-{}".format(length)
    return l

def get_merged_string(is_merged):
    m = ""
    if is_merged:
        m = ".merged-isoforms"
    return m

def get_micro_string(is_micro):
    m = ""
    if is_micro:
        m = ".micro-only"
    return m

def get_note_string(note=None):
    note_str = ""
    if (note is not None) and  (len(note) > 0):
        note_str = ".{}".format(note)
    return note_str

def get_offset_string(offset=None):
    o = ""
    if offset is not None:
        if isinstance(offset, (list, tuple)):
            o = "-".join(offset)
            o = ".offset-{}".format(o)
        else:
            o = ".offset-{}".format(offset)
    return o

def get_reweighting_iterations_string(reweighting_iterations=None):
    reweighting_iterations_str = ""
    if (reweighting_iterations is not None) and  (len(str(reweighting_iterations)) > 0):
        reweighting_iterations_str = ".rw-{}".format(reweighting_iterations)
    return reweighting_iterations_str


def get_smooth_string(is_smooth):
    s = ""
    if is_smooth:
        s = ".smooth"
    return s

def get_unique_string(is_unique):
    unique = ""
    if is_unique:
        unique = "-unique"
    return unique

def get_transcriptome_string(is_transcriptome):
    transcriptome = ""
    if is_transcriptome:
        transcriptome = ".transcriptome"
    return transcriptome


### b

def get_bed(base_path, name, is_merged=False):
    m = get_merged_string(is_merged)
    fn = '{}{}.bed.gz'.format(name, m)
    return os.path.join(base_path, fn)


def get_bf_rpkm_scatter_plot(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, fraction=None, is_smooth=False,
        reweighting_iterations=None,  note=None, image_type='pdf'):
    
    subfolder = os.path.join('orf-predictions', 'plots')
    s = get_riboseq_base(riboseq_base, name, subfolder, length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome,
        fraction=fraction, reweighting_iterations=reweighting_iterations, note=note)
    s = s + ".bf-rpkm.{}".format(image_type)
    return s


def get_bitseq_transcript_info(transcript_fasta):
    f = "{}.tr".format(transcript_fasta)
    return f

### c
def get_cds_bed(base_path, name, is_merged=False):
    m = get_merged_string(is_merged)
    fn = '{}.transcript-cds-coordinates{}.bed'.format(name, m)
    return os.path.join(base_path, 'transcript-index', fn)

def get_changepoint_image_file(base_path, note, condition, group, lookback, cp_threshold, image_type='pdf'):
    fn = "{}.{}.{}.lb-{}.cp-{}.{}".format(note, condition, group, lookback, cp_threshold, image_type)
    return os.path.join(base_path, 'plots', 'changepoints', fn)

### d
def get_diff_reg_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}-{}{}{}{}.diff-reg.{}'.format(condition_1, m, i, condition_2, m, i, n, image_type)
    return os.path.join(base_path, 'plots', 'diff-reg', fn)

def get_dominant_isoforms(base_path, note=None):
    n = get_note_string(note)
    fn = "dominant-isoforms{}.csv.gz".format(note)
    return os.path.join(base_path, fn)

### e
# used
def get_exons(base_path, name, note=None):
    note_str = get_note_string(note)
    fn = '{}.orfs-exons{}.bed.gz'.format(name, note_str)
    return os.path.join(base_path, 'transcript-index', fn)
 

### g

def get_gtf(base_path, name, is_merged=False):
    m = get_merged_string(is_merged)
    fn = '{}{}.gtf'.format(name, m)
    return os.path.join(base_path, fn)

### m

# used
def get_mean_and_var_image_file(base_path, condition, is_merged=False, is_isoforms=False, image_type='pdf', note=None):
    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}{}.mean-and-var.{}'.format(condition, m, i, n, image_type)
    return os.path.join(base_path, 'plots', 'mean-and-var', fn)


# used
def get_metagene_profiles(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, 
        is_transcriptome=False, is_merged=False, is_isoforms=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".metagene-profile.csv.gz"
    return s

# used
def get_metagene_profiles_bayes_factors(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_merged=False, is_isoforms=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".metagene-periodicity-bayes-factors.csv.gz"
    return s

# used
def get_default_models_base(project="rpbp_models"):
    import appdirs

    appname = "rpbp"
    appauthor = "dieterich-lab"
    models_base = appdirs.user_data_dir(appname, appauthor)
    models_base = os.path.join(models_base, project)
    return models_base


def get_models(models_base, model_type):
    path_ex = os.path.join(models_base, model_type, '*pkl')
    models = glob.glob(path_ex)
    return models

# used
def get_motif_analysis_base_folder(base, condition_1, condition_2, criterion, sequence_type):
    folder = "{}.{}.{}.{}".format(condition_1, condition_2, criterion, sequence_type)
    folder = os.path.join(base, 'motif-analysis', folder)
    return folder


def get_motif_analysis_folder(base, condition_1, condition_2, criterion, 
        sequence_type, fore_condition):

    folder = get_motif_analysis_base_folder(base, condition_1, condition_2, 
        criterion, sequence_type)

    subfolder = "{}.fore".format(fore_condition)
    return os.path.join(folder, subfolder)


def get_motif_analysis_results(base, condition_1, condition_2, criterion, 
        sequence_type, fore_condition):

    folder = get_motif_analysis_folder(base, condition_1, condition_2, 
        criterion, sequence_type, fore_condition)

    result_file = os.path.join(folder, "ame.txt")
    return result_file

def get_motif_sequences(base, condition_1, condition_2, criterion, direction, sequence_type):
    folder = get_motif_analysis_base_folder(base, condition_1, condition_2, 
        criterion, sequence_type)

    motif_sequences = "{}.{}-{}.{}.fa".format(condition_2, criterion, direction, sequence_type)
    motif_sequences = os.path.join(folder, motif_sequences)
    return motif_sequences

### o
# used
def get_orfs(base_path, name, note=None):
    note_str = get_note_string(note)
    fn = '{}.genomic-orfs{}.bed.gz'.format(name, note_str)
    return os.path.join(base_path, 'transcript-index', fn)


def get_orf_length_distribution_line_graph(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, fraction=None, is_smooth=False,
        reweighting_iterations=None,  note=None, is_grouped=False, is_chisq=False,
        image_type='pdf'):
    
    subfolder = os.path.join('orf-predictions', 'plots')
    s = get_riboseq_base(riboseq_base, name, subfolder, length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome,
        fraction=fraction, reweighting_iterations=reweighting_iterations, note=note, 
        is_grouped=is_grouped, is_chisq=is_chisq)
    s = s + ".orf-lengths.{}".format(image_type)
    return s


def get_orf_type_profile_base(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, fraction=None, is_smooth=False,
        reweighting_iterations=None,  note=None, is_grouped=False, is_chisq=False,
        subfolder='orf-predictions'):
    
    subfolder = os.path.join(subfolder, 'plots')
    s = get_riboseq_base(riboseq_base, name, subfolder, length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome,
        fraction=fraction, reweighting_iterations=reweighting_iterations, note=note, 
        is_grouped=is_grouped, is_chisq=is_chisq)
    return s



def get_orf_type_profile_image(base_path, orf_type, strand, image_type='eps'):
    fn = ".{}.{}.metagene-profiles.{}".format(orf_type, strand, image_type)
    return base_path + fn

def get_orf_types_pie_chart(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, fraction=None, is_smooth=False,
        reweighting_iterations=None,  note=None, is_grouped=False, is_chisq=False,
        is_filtered=False, image_type='pdf'):
    
    subfolder = os.path.join('orf-predictions', 'plots')
    s = get_riboseq_base(riboseq_base, name, subfolder, length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome,
        fraction=fraction, reweighting_iterations=reweighting_iterations, note=note, 
        is_grouped=is_grouped, is_chisq=is_chisq, is_filtered=is_filtered)
    s = s + ".orf-types.{}".format(image_type)
    return s


### p

# used
def get_peptide_coverage_line_graph(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, image_type='pdf'):
    
    subfolder = os.path.join('peptide-matches', 'plots')
    s = get_riboseq_base(riboseq_base, name, subfolder, length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        note=note)
    s = s + ".orf-peptide-coverage.{}".format(image_type)
    return s


# used
def get_periodic_offsets(riboseq_base, name, is_unique=False, is_cds_only=False, 
        is_transcriptome=False, is_merged=False, is_isoforms=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_merged=is_merged, is_isoforms=is_isoforms, note=note)
    s = s + ".periodic-offsets.csv.gz"
    return s

def get_preprocessing_report(base_path):
    fn = "preprocessing-report.tex"
    return os.path.join(base_path, 'preprocessing-report', fn)

### r
def get_raw_data_path(base_path):
    return os.path.join(base_path, 'raw-data')

def get_raw_data_fastqc_path(base_path):
    rdp = get_raw_data_path(base_path)
    return os.path.join(rdp, 'fastqc')


def get_raw_data_fastqc_data(base_path, filename):

    name = get_fastqc_name(filename)
    fastqc_folder = '{}_fastqc'.format(name)
    rdp = get_raw_data_fastqc_path(base_path)

    p = os.path.join(rdp, fastqc_folder, 'fastqc_data.txt')
    return p

### riboseq

# b

# used
def get_riboseq_base(riboseq_base, name, sub_folder, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_smooth=False, fraction=None, 
        reweighting_iterations=None, is_chisq=False, is_merged=False, is_isoforms=False, is_grouped=False, 
        is_filtered=False, note=None):
    
    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    o = get_offset_string(offset)
    transcriptome = get_transcriptome_string(is_transcriptome)
    chisq = get_chisq_string(is_chisq)
    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    s = get_smooth_string(is_smooth)
    f = get_fraction_string(fraction)
    r = get_reweighting_iterations_string(reweighting_iterations)
    g = get_grouped_string(is_grouped)
    fi = get_filtered_string(is_filtered)
    return os.path.join(riboseq_base, sub_folder, 
        '{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}'.format(
        name, n, transcriptome, m, i, unique, cds_only, l, o, s, f, r, chisq, g, fi))



# used

def get_riboseq_bam_base(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, 
        is_transcriptome=False, offset=None, is_merged=False, is_isoforms=False, note=None):
    
    bam_base = get_riboseq_base(riboseq_base, name, 'without-rrna-mapping', length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, offset=offset, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)
    return bam_base

# used
def get_riboseq_bam(riboseq_base, name, length=None, is_unique=False, is_cds_only=False, 
        is_transcriptome=False, offset=None, note=None, is_merged=False, is_isoforms=False):

    s = get_riboseq_bam_base(riboseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, offset=offset,
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)
    s = s + ".bam"
    return s

# used: get_all_read_filtering_counts
def get_riboseq_bam_fastqc_path(riboseq_data):
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'fastqc')

# used: get_all_read_filtering_counts
def get_riboseq_bam_fastqc_data(riboseq_data, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, is_chisq=False):

    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    n = get_note_string(note)
    c = get_chisq_string(is_chisq)
    name = '{}{}{}{}{}{}{}'.format(name, n, transcriptome, unique, cds_only, l, c)

    fastqc_folder = '{}_fastqc'.format(name)
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'fastqc', fastqc_folder, 'fastqc_data.txt')

# not used
def get_riboseq_bam_fastqc_read_lengths(riboseq_data, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, is_chisq=False):

    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    n = get_note_string(note)
    c = get_chisq_string(is_chisq)
    name = '{}{}{}{}{}{}{}'.format(name, n, transcriptome, unique, cds_only, l, c)

    fastqc_folder = '{}_fastqc'.format(name)
    return os.path.join(riboseq_data, 'without-rrna-mapping', 'fastqc', fastqc_folder, 
        'Images', 'sequence_length_distribution.png')


# used
def get_riboseq_bayes_factors(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_smooth=False, fraction=None, 
        reweighting_iterations=None, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
            is_smooth=is_smooth, fraction=fraction, reweighting_iterations=reweighting_iterations,
            note=note)
    
    s = s + ".bayes-factors.bed.gz"
    return s


def get_riboseq_bitseq(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, offset=None, is_merged=False, is_isoforms=False, note=None):

    unique = get_unique_string(is_unique)
    cds_only = get_cds_only_string(is_cds_only)
    l = get_length_string(length)
    o = get_offset_string(offset)
    m = get_merged_string(is_merged)   
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    transcriptome = get_transcriptome_string(is_transcriptome)
    
    return os.path.join(riboseq_base, 'transcript-abundance', 
        '{}{}{}{}{}{}{}{}{}.bitseq'.format(name, n, transcriptome, m, i, unique, cds_only, l, o))

def get_riboseq_bitseq_malphas(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, offset=None, is_merged=False, is_isoforms=False, note=None):

    s = get_riboseq_bitseq(riboseq_base, name, length=length, is_unique=is_unique, 
            is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".m_alphas"
    return s

def get_riboseq_bitseq_prob(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, offset=None, is_merged=False, is_isoforms=False, note=None):

    s = get_riboseq_bitseq(riboseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, offset=offset, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".prob"
    return s

def get_riboseq_bitseq_rpkm(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, offset=None, is_merged=False, is_isoforms=False, note=None):

    s = get_riboseq_bitseq(riboseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, offset=offset,
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".rpkm"
    return s

def get_riboseq_bitseq_rpkm_mean(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, offset=None, is_merged=False, is_isoforms=False, note=None):

    s = get_riboseq_bitseq(riboseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, offset=offset, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".rpkm.mean"
    return s

# f
def get_riboseq_fastq(riboseq_data, name):
    return os.path.join(riboseq_data, 'raw-data', '{}.fastq.gz'.format(name))

# m

def get_riboseq_metagene_profile_image(riboseq_base, name, image_type='eps', 
        length=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + "." + image_type
    return s

def get_metagene_profile_bayes_factor_image(riboseq_base, name, image_type='eps', 
        length=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', length=length, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".bayes-factors." + image_type
    return s



# p

# used
def get_riboseq_peptide_matches(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'peptide-matches', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_chisq=is_chisq, note=note)
    s = s + ".peptide-matches.csv.gz"
    return s


# used
def get_riboseq_predicted_orfs(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, is_smooth=False, fraction=None, 
        reweighting_iterations=None, is_chisq=False, is_filtered=False):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_smooth=is_smooth, fraction=fraction, reweighting_iterations=reweighting_iterations, 
        is_chisq=is_chisq, is_filtered=is_filtered, note=note)
    s = s + ".predicted-orfs.bed.gz"
    return s

# used
def get_riboseq_predicted_orfs_dna(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, is_smooth=False, fraction=None, 
        reweighting_iterations=None, is_filtered=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_smooth=is_smooth, fraction=fraction, reweighting_iterations=reweighting_iterations, 
        is_chisq=is_chisq, is_filtered=is_filtered, note=note)
    s = s + ".predicted-orfs.dna.fa"
    return s

# used
def get_riboseq_predicted_orfs_protein(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_chisq=False, is_smooth=False, fraction=None, 
        reweighting_iterations=None,is_filtered=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome,
        is_smooth=is_smooth, fraction=fraction, reweighting_iterations=reweighting_iterations,  
        is_chisq=is_chisq, is_filtered=is_filtered, note=note)
    s = s + ".predicted-orfs.protein.fa"
    return s

def get_riboseq_predicted_orf_mackowiak_overlap_image(riboseq_base, name, image_type='eps', 
        length=None, offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, 
        is_chisq=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_chisq=is_chisq, note=note)
    s = s + ".predicted-orf-mackowiak-overlap.{}".format(image_type)
    return s

def get_riboseq_predicted_orf_peptide_coverage(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-peptide-coverage.csv.gz"
    return s


def get_riboseq_predicted_orf_peptide_coverage_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-peptide-coverage.{}".format(image_type)
    return s

def get_riboseq_predicted_orf_qti_seq_overlap_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-qti-seq-overlap.{}".format(image_type)
    return s


def get_riboseq_predicted_orf_spikes(riboseq_base, name, length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".spike-codons.csv.gz"
    return s
    
def get_riboseq_predicted_orf_spikes_bed(riboseq_base, name, length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".spikes.bed.gz"
    return s

def get_riboseq_predicted_orf_spikes_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".spike-codons.{}".format(image_type)
    return s



def get_riboseq_predicted_orf_type_overlap_image(riboseq_base, name, image_type='eps', length=None, 
        offset=None, is_unique=False, is_cds_only=False, is_transcriptome=False, note=None):
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', length=length, offset=offset, 
        is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)
    s = s + ".predicted-orf-type-overlap.{}".format(image_type)
    return s


# used
def get_riboseq_profiles(riboseq_base, name, length=None, offset=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_smooth=False, fraction=None, 
        reweighting_iterations=None, note=None):

    s = get_riboseq_base(riboseq_base, name, 'orf-profiles', length=length, offset=offset, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
            is_smooth=is_smooth, fraction=fraction, 
            reweighting_iterations=reweighting_iterations, note=note)

    s = s + ".profiles.mtx.gz"
    return s

# r

def get_riboseq_read_filtering_counts(riboseq_base, note=None):
    note_str = get_note_string(note)
    fn = "read-filtering-counts{}.csv.gz".format(note_str)
    s = os.path.join(riboseq_base, fn)
    return s

def get_riboseq_read_filtering_counts_image(riboseq_base, note="", image_type="eps"):
    note_str = get_note_string(note)
    fn = "read-filtering-counts{}.{}".format(note_str, image_type)
    s = os.path.join(riboseq_base, fn)
    return s

def get_riboseq_read_length_distribution(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, offset=None, note=None):

    s = get_riboseq_base(riboseq_base, name, 'without-rrna-mapping', length=length, offset=offset, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)

    s = s + ".length-distribution.csv.gz"
    return s


def get_riboseq_read_length_distribution_image(riboseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, note=None, offset=None, image_type='eps'):

    subfolder = os.path.join('without-rrna-mapping', 'plots')

    s = get_riboseq_base(riboseq_base, name, subfolder, length=length, offset=offset, 
            is_unique=is_unique, is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, note=note)

    s = s + ".length-distribution.{}".format(image_type)
    return s


### rna

# b
def get_rnaseq_bam_base(rnaseq_base, name, length=None, is_unique=False, is_cds_only=False, 
        is_transcriptome=False, is_merged=False, is_isoforms=False, note=None):

    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)

    return os.path.join(rnaseq_base, 'mapping', '{}{}{}{}{}{}{}{}'.format(name, n, transcriptome, 
        m, i, unique, cds_only, l))


def get_rnaseq_bam(rnaseq_base, name, length=None, is_unique=False, is_cds_only=False, 
        is_transcriptome=False, is_merged=False, is_isoforms=False, note=None):
    
    s = get_rnaseq_bam_base(rnaseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)
    s = s + ".bam"
    return s

def get_rnaseq_bam_path(base_path):
    return os.path.join(base_path, 'mapping',)

def get_rnaseq_bitseq(rnaseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_merged=False, 
        is_isoforms=False, note=None):
    
    unique = get_unique_string(is_unique)
    cds_only = get_cds_only_string(is_cds_only)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)

    n = get_note_string(note)

    return os.path.join(rnaseq_base, 'transcript-abundance', 
        '{}-rna{}{}{}{}{}{}{}.bitseq'.format(name, n, transcriptome, m, i, unique, cds_only, l))

def get_rnaseq_bitseq_malphas(rnaseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_merged=False, 
        is_isoforms=False, note=None):

    s = get_rnaseq_bitseq(rnaseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".m_alphas"
    return s

def get_rnaseq_bitseq_prob(rnaseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_merged=False, 
        is_isoforms=False, note=None):

    s = get_rnaseq_bitseq(rnaseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".prob"
    return s

def get_rnaseq_bitseq_rpkm(rnaseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_merged=False, 
        is_isoforms=False, note=None):

    s = get_rnaseq_bitseq(rnaseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".rpkm"
    return s


def get_rnaseq_bitseq_rpkm_mean(rnaseq_base, name, length=None, is_unique=False, 
        is_cds_only=False, is_transcriptome=False, is_merged=False, 
        is_isoforms=False, note=None):

    s = get_rnaseq_bitseq(rnaseq_base, name, length=length, is_unique=is_unique, 
        is_cds_only=is_cds_only, is_transcriptome=is_transcriptome, 
        is_merged=is_merged, is_isoforms=is_isoforms, note=note)

    s = s + ".rpkm.mean"
    return s


# used
def get_rpkm_image_file(base_path, condition, is_merged=False, is_isoforms=False, image_type='pdf', note=None):
    m = get_merged_string(is_merged)
    n = get_note_string(note)
    fn = '{}{}{}.rpkm.{}'.format(condition, m, n, image_type)
    return os.path.join(base_path, 'plots', 'rpkm', fn)

# used
def get_rpkm_fold_change_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', note=None):

    m = get_merged_string(is_merged)
    n = get_note_string(note)
    fn = '{}{}-{}{}{}.rpkm-fc.{}'.format(condition_1, m, condition_2, m, n, image_type)
    return os.path.join(base_path, 'plots', 'rpkm-fc', fn)

# used
def get_rpkm_vs_rpkm_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', field='ribo', note=None):

    m = get_merged_string(is_merged)
    n = get_note_string(note)
    fn = '{}{}-{}{}{}.{}-rpkm-vs-rpkm.{}'.format(condition_1, m, condition_2, m, n, field, image_type)
    return os.path.join(base_path, 'plots', 'rpkm-vs-rpkm', fn)

### s

# used
def get_star_index(base_path, name, is_merged=False):
    m = get_merged_string(is_merged)
    fn = '{}{}'.format(name, m)
    return os.path.join(base_path, 'STAR', fn)


### t
def get_te_kl(base_path, name, is_merged=False, is_isoforms=False, note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}{}.te-kl.csv.gz'.format(name, m, i, n)
    return os.path.join(base_path, fn)


def get_te_kl_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}-{}{}{}{}.te-kl.{}'.format(condition_1, m, i, condition_2, m, i, n, image_type)
    return os.path.join(base_path, 'plots', 'te-kl', fn)


def get_te_ma_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}-{}{}{}{}.te-ma.{}'.format(condition_1, m, i, condition_2, m, i, n, image_type)
    return os.path.join(base_path, 'plots', 'te-ma', fn)

def get_te_pvalues(base_path, name, is_merged=False, is_isoforms=False, note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}{}.te-pvalues.csv.gz'.format(name, m, i, n)
    return os.path.join(base_path, fn)


def get_te_rpkm_fold_change_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}-{}{}{}{}.te-rpkm-fc.{}'.format(condition_1, m, i, condition_2, m, i, n, image_type)
    return os.path.join(base_path, 'plots', 'te-rpkm-fc', fn)


def get_transcript_fasta(base_path, name, is_merged=False):
    m = get_merged_string(is_merged)
    fn = '{}.transcripts{}.fa'.format(name, m)
    return os.path.join(base_path, 'transcript-index', fn)

def get_translational_efficiency(base_path, condition, is_merged=False, is_isoforms=False, note=None):
    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}{}.translational-efficiency.csv.gz'.format(condition, m, i, n)
    return os.path.join(base_path, fn)

### w

# used
def get_without_adapters_base(base_path, name, note=None):
    n = get_note_string(note)
    base = "{}{}".format(name, n)
    return os.path.join(base_path, 'without-adapters', base)

# used
def get_without_adapters_fastq(base_path, name, note=None):
    base = get_without_adapters_base(base_path, name, note=note)
    fastq = "{}.fastq.gz".format(base)
    return fastq

def get_without_adapters_fastqc(base_path):
    return os.path.join(base_path, 'without-adapters', 'fastqc')

def get_without_adapters_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = '{}{}_fastqc'.format(name, n)
    return os.path.join(base_path, 'without-adapters', 'fastqc', fastqc_folder, 'fastqc_data.txt')


# used
def get_with_rrna_fastq(base_path, name, note=None):
    n = get_note_string(note)
    name = "{}{}".format(name, n)
    return os.path.join(base_path, 'with-rrna', '{}.fastq.gz'.format(name))

# used
def get_with_rrna_fastqc(base_path):
    return os.path.join(base_path, 'with-rrna', 'fastqc')

def get_with_rrna_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = '{}{}_fastqc'.format(name, n)
    return os.path.join(base_path, 'with-rrna', 'fastqc', fastqc_folder, 'fastqc_data.txt')

# used
def get_without_rrna_fastq(base_path, name, note=None):
    n = get_note_string(note)
    name = "{}{}".format(name, n)
    return os.path.join(base_path, 'without-rrna', '{}.fastq.gz'.format(name))

def get_without_rrna_fastqc(base_path):
    return os.path.join(base_path, 'without-rrna', 'fastqc')

def get_without_rrna_fastqc_data(base_path, name, note=None):
    n = get_note_string(note)
    fastqc_folder = '{}{}_fastqc'.format(name, n)
    return os.path.join(base_path, 'without-rrna', 'fastqc', fastqc_folder, 'fastqc_data.txt')

