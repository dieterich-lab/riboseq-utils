import glob
import os

import misc.utils as utils

### parameterized names

def get_annotated_string(is_annotated):
    annotated = ""
    if is_annotated:
        annotated = ".annotated"
    return annotated

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

def get_de_novo_string(is_de_novo):
    de_novo = ""
    if is_de_novo:
        de_novo = ".de-novo"
    return de_novo

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

def get_isoform_strategy_string(isoform_strategy):
    s = ""
    if (isoform_strategy is not None) and (len(isoform_strategy) > 0):
        s = ".{}".format(isoform_strategy)

    return s


def get_length_string(length=None):
    l = ""
    if length is not None:
        if isinstance(length, (list, tuple)):
            l = "-".join(str(l) for l in length)
            l = ".length-{}".format(l)
        else:
            l = ".length-{}".format(str(length))
    return l

def get_merged_string(is_merged):
    m = ""
    if is_merged:
        m = ".merged-isoforms"
    return m

def get_star_input_string(is_star_input):
    s = ""
    if is_star_input:
        s = ".star-input"
    return s

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
            o = "-".join(str(o) for o in offset)
            o = ".offset-{}".format(o)
        else:
            o = ".offset-{}".format(str(offset))
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

def get_zscore_string(is_zscore):
    s = ""
    if is_zscore:
        s = ".zscore"
    return s

def get_transcriptome_string(is_transcriptome):
    transcriptome = ""
    if is_transcriptome:
        transcriptome = ".transcriptome"
    return transcriptome

def get_unique_string(is_unique):
    unique = ""
    if is_unique:
        unique = "-unique"
    return unique

### a

def get_abundances(base_path, isoform_strategy=None, note=None):
    n = get_note_string(note)
    i = get_isoform_strategy_string(isoform_strategy)
    fn = "abundances{}{}.csv.gz".format(i, n)
    return os.path.join(base_path, fn)

### b

def get_b_tea_differential_analysis_report(base_path, note):
    n = get_note_string(note)
    fn = "differential-analysis-report{}.tex".format(n)
    return os.path.join(base_path, fn)


def get_bed(
        base_path, 
        name, 
        is_merged=False, 
        is_annotated=False, 
        is_de_novo=False, 
        is_cds_only=False
    ):

    m = get_merged_string(is_merged)
    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    cds = get_cds_only_string(is_cds_only)
    fn = '{}{}{}{}{}.bed.gz'.format(name, m, c, d, cds)
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

def get_differential_gene_xlsx(
        base_path,
        condition_1,
        condition_2, 
        isoform_strategy=None,
        is_zscore=False,
        note=None):


    fn = [
        condition_1,
        "-",
        condition_2,
        get_isoform_strategy_string(isoform_strategy),
        get_zscore_string(is_zscore),
        get_note_string(note),
        ".diff-genes.xlsx"
    ]
    fn = ''.join(fn)
    return os.path.join(base_path, 'diff-genes', fn)


def get_diff_reg_image_file(
        base_path,
        condition_1,
        condition_2, 
        isoform_strategy=None,
        is_zscore=False,
        image_type='pdf',
        note=None):

    
    fn = [
        condition_1,
        "-",
        condition_2,
        get_isoform_strategy_string(isoform_strategy),
        get_zscore_string(is_zscore),
        get_note_string(note),
        ".diff-reg.",
        image_type
    ]
    fn = ''.join(fn)
    return os.path.join(base_path, 'plots', 'diff-reg', fn)

### e
# used
def get_exons(
        base_path, 
        name, 
        is_annotated=False, 
        is_de_novo=False,
        is_orf=False,
        note=None
    ):

    note_str = get_note_string(note)
    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    orf_subscr = 'orfs' if is_orf else 'exons'
    fn = '{}.{}{}{}{}.bed.gz'.format(name, orf_subscr, c, d, note_str)
    return os.path.join(base_path, 'transcript-index', fn)
 

### g

def get_gtf(
        base_path, 
        name, 
        is_de_novo=False, 
        is_annotated=False, 
        is_merged=False, 
        is_cds_only=False,
        is_gff3=False,
        is_star_input=False
    ):

    c = get_annotated_string(is_annotated)
    m = get_merged_string(is_merged)
    cds = get_cds_only_string(is_cds_only)
    d = get_de_novo_string(is_de_novo)
    s = get_star_input_string(is_star_input)
    ext = 'gff' if is_gff3 else 'gtf'
    fn = '{}{}{}{}{}{}.{}'.format(name, m, c, d, cds, s, ext)
    return os.path.join(base_path, fn)

### m

# used
def get_mean_and_var_image_file(
        base_path, 
        condition,
        isoform_strategy=None,
        image_type='pdf',
        note=None):

    fn = [
        condition,
        get_isoform_strategy_string(isoform_strategy),
        get_note_string(note),
        ".mean-and-var.",
        image_type
    ]
    fn = ''.join(fn)
    return os.path.join(base_path, 'plots', 'mean-and-var', fn)

# used
def get_metagene_profiles(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', **kwargs)
    s = s + ".metagene-profile.csv.gz"
    return s

# used
def get_metagene_profiles_bayes_factors(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', **kwargs)
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
    import shlex

    path_ex = os.path.join(models_base, model_type, '*pkl')
    models = glob.glob(path_ex)
    models = [shlex.quote(m) for m in models]
    return models

# used
def get_motif_analysis_base_folder(
        base,
        condition_1,
        condition_2,
        field,
        region,
        isoform_strategy):

    folder = [
        condition_1,
        condition_2,
        field,
        region,
        isoform_strategy
    ]

    folder = '.'.join(folder)
    folder = os.path.join(base, 'motif-analysis', folder)
    return folder


def get_motif_analysis_folder(fore_condition, **kwargs):
    folder = get_motif_analysis_base_folder(**kwargs)
    subfolder = "{}.fore".format(fore_condition)
    return os.path.join(folder, subfolder)


def get_motif_analysis_results(**kwargs):
    folder = get_motif_analysis_folder(**kwargs)
    result_file = os.path.join(folder, "ame.txt")
    return result_file

def get_motif_fimo_folder(motifs, direction, **kwargs):
    
    base_folder = get_motif_analysis_base_folder(**kwargs)
    motifs_name = utils.get_basename(motifs)

    fimo_folder = [
        kwargs['condition_2'],
        direction,
        motifs_name
    ]

    fimo_folder = '.'.join(fimo_folder)
    fimo_folder = os.path.join(base_folder, fimo_folder)
    return fimo_folder


def get_motif_fimo_results(**kwargs):   
    fimo_folder = get_motif_fimo_folder(**kwargs)
    result_file = os.path.join(fimo_folder, "fimo.txt")
    return result_file

def get_motif_heatmap(
        base,
        condition_1,
        condition_2,
        isoform_strategy,
        image_type="pdf"):
    
    folder = [
        condition_1,
        condition_2,
        isoform_strategy
    ]

    folder = '.'.join(folder)
    folder = os.path.join(base, 'motif-analysis', folder)

    fn = "overrepresented-motifs-heatmap.{}".format(image_type)
    return os.path.join(folder, fn)

def get_motif_report(base, condition_1, condition_2, isoform_strategy):
    
    folder = [
        condition_1,
        condition_2,
        isoform_strategy
    ]

    folder = '.'.join(folder)
    folder = os.path.join(base, 'motif-analysis', folder)

    fn = "overrepresented-motifs.tex"
    return os.path.join(folder, fn)

def get_motif_results(base, isoform_strategy):
    i = get_isoform_strategy_string(isoform_strategy)
    motif_results = [
        "overrepresented-motifs",
        i,
        ".csv.gz"
    ]
    motif_results = ''.join(motif_results)
    motif_results = os.path.join(base, 'motif-analysis', motif_results)
    return motif_results

    
def get_motif_sequences(direction, **kwargs):
       
    folder = get_motif_analysis_base_folder(**kwargs)

    motif_sequences = [
        kwargs['condition_2'],
        direction
    ]
    motif_sequences = '.'.join(motif_sequences)
    motif_sequences = motif_sequences + ".fa"
    motif_sequences = os.path.join(folder, motif_sequences)
    return motif_sequences

### o
# used
def get_orfs(base_path, name, is_annotated=False, is_de_novo=False, note=None):
    note_str = get_note_string(note)
    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    fn = '{}.genomic-orfs{}{}{}.bed.gz'.format(name, c, d, note_str)
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

def get_orf_types_pie_chart(
    riboseq_base, 
    name, 
    image_type='pdf',
    **kwargs):
    
    subfolder = os.path.join('orf-predictions', 'plots')
    s = get_riboseq_base(
        riboseq_base, 
        name, 
        subfolder,
        **kwargs
    )
    s = s + ".orf-types-pie.{}".format(image_type)
    return s


def get_orf_types_bar_chart(
    riboseq_base, 
    name, 
    image_type='pdf',
    **kwargs):
    
    subfolder = os.path.join('orf-predictions', 'plots')
    s = get_riboseq_base(
        riboseq_base, 
        name, 
        subfolder,
        **kwargs
    )
    s = s + ".orf-types-bar.{}".format(image_type)
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
def get_periodic_offsets(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', **kwargs)
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

### ribodiff
def get_ribodiff_base(base, condition_1, condition_2, is_unique=False):
    unique = get_unique_string(is_unique)
    subfolder = "{}_{}{}".format(condition_1, condition_2, unique)
    return os.path.join(base, 'ribodiff', subfolder, "{}_{}".format(condition_1, condition_2))

def get_ribodiff_design(base, condition_1, condition_2, is_unique=False):
    base = get_ribodiff_base(base, condition_1, condition_2, is_unique=is_unique)
    f = base + ".experimental_design.csv"
    return f

def get_ribodiff_gene_count_table(base, condition_1, condition_2, is_unique=False):
    base = get_ribodiff_base(base, condition_1, condition_2, is_unique=is_unique)
    f = base + ".gene_count_table.tab"
    return f

def get_ribodiff_results(base, condition_1, condition_2, is_unique=False):
    base = get_ribodiff_base(base, condition_1, condition_2, is_unique=is_unique)
    f = base + ".results.tab"
    return f

### call either rna or riboseq bam file

def get_seq_bam(seq, base, name, **kwargs):
    s = ''
    if seq == 'rna':
        s = get_rnaseq_bam(base, name, **kwargs)
    if seq == 'ribo':
        s = get_riboseq_bam(base, name, **kwargs)
    return s


### riboseq

# b

# used
def get_riboseq_base(
        riboseq_base,
        name,
        sub_folder,
        length=None,
        offset=None,
        is_unique=False, 
        is_cds_only=False,
        is_transcriptome=False,
        is_smooth=False,
        fraction=None, 
        reweighting_iterations=None,
        is_chisq=False,
        isoform_strategy=None,
        is_grouped=False, 
        is_filtered=False,
        note=None):
    
    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    o = get_offset_string(offset)
    transcriptome = get_transcriptome_string(is_transcriptome)
    chisq = get_chisq_string(is_chisq)
    i = get_isoform_strategy_string(isoform_strategy)
    n = get_note_string(note)
    s = get_smooth_string(is_smooth)
    f = get_fraction_string(fraction)
    r = get_reweighting_iterations_string(reweighting_iterations)
    g = get_grouped_string(is_grouped)
    fi = get_filtered_string(is_filtered)

    fn = ''.join([
        name,
        n,
        transcriptome,
        i,
        unique,
        cds_only,
        l,
        o,
        s,
        f,
        r,
        chisq,
        g,
        fi
    ])

    return os.path.join(riboseq_base, sub_folder, fn)
        



# used

def get_riboseq_bam_base(riboseq_base, name, **kwargs):
    
    bam_base = get_riboseq_base(
        riboseq_base,
        name,
        'without-rrna-mapping',
        **kwargs
    )

    return bam_base

# used
def get_riboseq_bam(riboseq_base, name, **kwargs):

    s = get_riboseq_bam_base(riboseq_base, name, **kwargs)
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
def get_riboseq_bayes_factors(riboseq_base, name, **kwargs):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    
    s = s + ".bayes-factors.bed.gz"
    return s


def get_riboseq_bitseq(
        riboseq_base,
        name,
        length=None,
        is_unique=False, 
        is_cds_only=False,
        is_transcriptome=False,
        offset=None,
        isoform_strategy=None,
        note=None):

    unique = get_unique_string(is_unique)
    cds_only = get_cds_only_string(is_cds_only)
    l = get_length_string(length)
    o = get_offset_string(offset)
    i = get_isoform_strategy_string(isoform_strategy)
    n = get_note_string(note)
    transcriptome = get_transcriptome_string(is_transcriptome)
    
    fn = ''.join([
        name,
        n,
        transcriptome,
        i,
        unique,
        cds_only,
        l,
        o,
        '.bitseq'
    ])
    return os.path.join(riboseq_base, 'transcript-abundance', fn)

def get_riboseq_bitseq_malphas(riboseq_base, name, **kwargs):

    s = get_riboseq_bitseq(riboseq_base, name, **kwargs)
    s = s + ".m_alphas"
    return s

def get_riboseq_bitseq_prob(riboseq_base, name, **kwargs):

    s = get_riboseq_bitseq(riboseq_base, name, **kwargs)
    s = s + ".prob"
    return s

def get_riboseq_bitseq_rpkm(riboseq_base, name, **kwargs):

    s = get_riboseq_bitseq(riboseq_base, name, **kwargs)
    s = s + ".rpkm"
    return s

def get_riboseq_bitseq_rpkm_mean(riboseq_base, name, **kwargs):

    s = get_riboseq_bitseq(riboseq_base, name, **kwargs)
    s = s + ".rpkm.mean"
    return s

# c
def get_riboseq_cell_type_protein(riboseq_base, name, **kwargs):
    
    s = get_riboseq_base(riboseq_base, name, 'cell-types', **kwargs)
    s = s + ".predicted-orfs.protein.fa"
    return s


# f
def get_riboseq_fastq(riboseq_data, name):
    return os.path.join(riboseq_data, 'raw-data', '{}.fastq.gz'.format(name))

# m

def get_metagene_profile_image(base, name, image_type='eps', **kwargs):

    s = get_riboseq_base(base, name, 'metagene-profiles', **kwargs)
    s = s + "." + image_type
    return s

def get_metagene_profile_bayes_factor_image(
        riboseq_base,
        name,
        image_type='eps',
        **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'metagene-profiles', **kwargs)
    s = s + ".bayes-factors." + image_type
    return s



# p

# used
def get_riboseq_peptide_matches(riboseq_base, name, peptide_name, **kwargs):
    
    n = "{}-{}".format(name, peptide_name)

    s = get_riboseq_base(riboseq_base, n, 'peptide-matches', **kwargs)
    s = s + ".peptide-matches.csv.gz"
    return s


# used
def get_riboseq_predicted_orfs(riboseq_base, name, **kwargs):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    s = s + ".predicted-orfs.bed.gz"
    return s

# used
def get_riboseq_predicted_orfs_dna(riboseq_base, name, **kwargs):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    s = s + ".predicted-orfs.dna.fa"
    return s

# used
def get_riboseq_predicted_orfs_protein(riboseq_base, name, **kwargs):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    s = s + ".predicted-orfs.protein.fa"
    return s

def get_riboseq_predicted_orf_peptide_coverage(riboseq_base, name, **kwargs):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    s = s + ".predicted-orf-peptide-coverage.csv.gz"
    return s


def get_riboseq_predicted_orf_peptide_coverage_image(
        riboseq_base,
        name,
        image_type='eps',
        **kwargs):
    
    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    s = s + ".predicted-orf-peptide-coverage.{}".format(image_type)
    return s

def get_riboseq_predicted_orf_qti_seq_overlap_image(
        riboseq_base,
        name,
        image_type='eps',
        **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    s = s + ".predicted-orf-qti-seq-overlap.{}".format(image_type)
    return s

def get_riboseq_predicted_orf_type_overlap_image(
        riboseq_base,
        name,
        image_type='eps',
        **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'orf-predictions', **kwargs)
    s = s + ".predicted-orf-type-overlap.{}".format(image_type)
    return s


# used
def get_riboseq_profiles(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'orf-profiles', **kwargs)
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

def get_riboseq_read_length_distribution(riboseq_base, name, **kwargs):

    s = get_riboseq_base(riboseq_base, name, 'without-rrna-mapping', **kwargs)
    s = s + ".length-distribution.csv.gz"
    return s


def get_riboseq_read_length_distribution_image(riboseq_base, name, image_type='eps', **kwargs):

    subfolder = os.path.join('without-rrna-mapping', 'plots')

    s = get_riboseq_base(riboseq_base, name, subfolder, **kwargs)
    s = s + ".length-distribution.{}".format(image_type)
    return s


### rna

# b
def get_rnaseq_bam_base(
        rnaseq_base,
        name,
        length=None,
        is_unique=False,
        is_cds_only=False,
        is_transcriptome=False,
        isoform_strategy=None,
        note=None):

    cds_only = get_cds_only_string(is_cds_only)
    unique = get_unique_string(is_unique)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    i = get_isoform_strategy_string(isoform_strategy)
    n = get_note_string(note)

    bam_base = '{}{}{}{}{}{}{}'.format(
        name, 
        n, 
        transcriptome, 
        i, 
        unique, 
        cds_only, 
        l
    )

    rnaseq_bam_path = get_rnaseq_bam_path(rnaseq_base)
    bam_base = os.path.join(rnaseq_bam_path, bam_base)
    return bam_base


def get_rnaseq_bam(rnaseq_base, name, **kwargs):
    
    s = get_rnaseq_bam_base(rnaseq_base, name, **kwargs)
    s = s + ".bam"
    return s

def get_rnaseq_bam_path(base_path):
    return os.path.join(base_path, 'mapping',)

def get_rnaseq_bitseq(
        rnaseq_base,
        name,
        length=None,
        is_unique=False, 
        is_cds_only=False,
        is_transcriptome=False,
        isoform_strategy=None,
        note=None):
    
    unique = get_unique_string(is_unique)
    cds_only = get_cds_only_string(is_cds_only)
    l = get_length_string(length)
    transcriptome = get_transcriptome_string(is_transcriptome)
    i = get_isoform_strategy_string(isoform_strategy)
    n = get_note_string(note)

    fn = '{}-rna{}{}{}{}{}{}.bitseq'.format(
        name,
        n,
        transcriptome,
        i,
        unique,
        cds_only,
        l
    )

    return os.path.join(rnaseq_base, 'transcript-abundance', fn)

def get_rnaseq_bitseq_malphas(rnaseq_base, name, **kwargs):

    s = get_rnaseq_bitseq(rnaseq_base, name, **kwargs)
    s = s + ".m_alphas"
    return s

def get_rnaseq_bitseq_prob(rnaseq_base, name, **kwargs):

    s = get_rnaseq_bitseq(rnaseq_base, name, **kwargs)
    s = s + ".prob"
    return s

def get_rnaseq_bitseq_rpkm(rnaseq_base, name, **kwargs):

    s = get_rnaseq_bitseq(rnaseq_base, name, **kwargs)
    s = s + ".rpkm"
    return s


def get_rnaseq_bitseq_rpkm_mean(rnaseq_base, name, **kwargs):

    s = get_rnaseq_bitseq(rnaseq_base, name, **kwargs)
    s = s + ".rpkm.mean"
    return s

# used
def get_rpbp_prediction_report(base_path, note):
    n = get_note_string(note)
    fn = "prediction-report{}.tex".format(n)
    return os.path.join(base_path, fn)

# used
def get_rpkm_image_file(
        base_path,
        condition,
        isoform_strategy=None,
        image_type='pdf',
        note=None):

    fn = [
        condition,
        get_isoform_strategy_string(isoform_strategy),
        get_note_string(note),
        ".rpkm.",
        image_type
    ]
    fn = ''.join(fn)

    return os.path.join(base_path, 'plots', 'rpkm', fn)

# used
def get_rpkm_fold_change_image_file(
        base_path, 
        condition_1, 
        condition_2, 
        isoform_strategy=None,
        is_filtered=False,
        is_zscore=False,
        image_type='pdf',
        note=None):
    
    fn = [
        condition_1,
        "-",
        condition_2,
        get_isoform_strategy_string(isoform_strategy),
        get_zscore_string(is_zscore),
        get_filtered_string(is_filtered),
        get_note_string(note),
        ".rpkm-fc.",
        image_type
    ]
    fn = ''.join(fn)
    return os.path.join(base_path, 'plots', 'rpkm-fc', fn)

# used
def get_rpkm_vs_rpkm_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', field='ribo', note=None):

    m = get_merged_string(is_merged)
    n = get_note_string(note)
    fn = '{}{}-{}{}{}.{}-rpkm-vs-rpkm.{}'.format(condition_1, m, condition_2, m, n, field, image_type)
    return os.path.join(base_path, 'plots', 'rpkm-vs-rpkm', fn)

def get_rpkm_te_comparison_image_file(
        base_path,
        condition_1,
        condition_2,
        isoform_strategy=None,
        is_filtered=False,
        is_zscore=False,
        image_type='pdf',
        note=None):

    fn = [
        condition_1,
        "-",
        condition_2,
        get_isoform_strategy_string(isoform_strategy),
        get_zscore_string(is_zscore),
        get_filtered_string(is_filtered),
        get_note_string(note),
        ".rpkm-te-comparison.",
        image_type
    ]
    fn = ''.join(fn)

    return os.path.join(base_path, 'plots', 'rpkm-te-comparison', fn)


### s

# used
def get_sample_embedding_file(
        base_path,
        isoform_strategy=None,
        image_type='pdf',
        note=None):
     
    fn = [
        "sample-embedding",
        get_isoform_strategy_string(isoform_strategy),
        get_note_string(note),
        ".",
        image_type
    ]

    fn = ''.join(fn)
    return os.path.join(base_path, 'plots', fn)

def get_star_index(base_path, name, is_merged=False):
    m = get_merged_string(is_merged)
    fn = '{}{}'.format(name, m)
    return os.path.join(base_path, 'STAR', fn)


### t
def get_te_kl(base_path, name, isoform_strategy=None, note=None):

    i = get_isoform_strategy_string(isoform_strategy)
    n = get_note_string(note)
    fn = ''.join([
        name,
        i,
        n,
        '.te-kl.csv.gz'
    ])
    return os.path.join(base_path, fn)


def get_te_kl_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}-{}{}{}{}.te-kl.{}'.format(condition_1, m, i, condition_2, m, i, n, image_type)
    return os.path.join(base_path, 'plots', 'te-kl', fn)


def get_ma_image_file(
        base_path,
        condition_1,
        condition_2, 
        isoform_strategy=None,
        is_filtered=False,
        is_zscore=False,
        image_type='pdf',
        note=None):

    
    fn = [
        condition_1,
        "-",
        condition_2,
        get_isoform_strategy_string(isoform_strategy),
        get_zscore_string(is_zscore),
        get_filtered_string(is_filtered),
        get_note_string(note),
        ".ma-plots.",
        image_type
    ]
    fn = ''.join(fn)

    return os.path.join(base_path, 'plots', 'ma-plots', fn)

def get_te_pvalues(
        base_path, 
        name, 
        isoform_strategy=None,
        is_zscore=False,
        note=None):

    i = get_isoform_strategy_string(isoform_strategy)
    z = get_zscore_string(is_zscore)
    n = get_note_string(note)

    fn = ''.join([
        name,
        i,
        n,
        ".te-pvalues",
        z,
        ".csv.gz"
    ])
    return os.path.join(base_path, fn)

def get_te_rpkm_fold_change_image_file(base_path, condition_1, condition_2, 
        is_merged=False, is_isoforms=False, image_type='pdf', note=None):

    m = get_merged_string(is_merged)
    i = get_isoforms_string(is_isoforms)
    n = get_note_string(note)
    fn = '{}{}{}-{}{}{}{}.te-rpkm-fc.{}'.format(condition_1, m, i, condition_2, m, i, n, image_type)
    return os.path.join(base_path, 'plots', 'te-rpkm-fc', fn)


def get_transcript_fasta(
        base_path,
        name,
        is_merged=False,
        is_annotated=False,
        is_de_novo=False,
        is_cds_only=False
    ):

    m = get_merged_string(is_merged)
    c = get_annotated_string(is_annotated)
    d = get_de_novo_string(is_de_novo)
    cds = get_cds_only_string(is_cds_only)
    fn = '{}.transcripts{}{}{}{}.fa'.format(name, m, c, d, cds)
    return os.path.join(base_path, 'transcript-index', fn)

def get_translational_efficiency(
        base_path,
        condition,
        isoform_strategy=None,
        note=None):

    i = get_isoform_strategy_string(isoform_strategy)
    n = get_note_string(note)

    fn = ''.join([
        condition,
        i,
        n,
        ".translational-efficiency.csv.gz"
    ])
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

