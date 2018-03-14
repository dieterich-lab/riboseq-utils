######################################################################
# This file is an ad hoc script used to simulate data for the b-tea
# pipeline. TE can be estimated using a simulation approach based on 
# modifying real data to introduce differences in RNA/RFP data. 
# This script uses samtools to subsample existing RNA/RFP BAM files
# into 3 groups:
#
#   (i) RNA + RFP: concordant differences in expression that
#   should be detectable in both RNA and RFP; 
#   (ii) RFP only: translationally induced differences which 
#   are exclusive to RFP;
#   (iii) RNA only: "buffered" genes that do not exhibit variation
#   on the translational level, but that should be detected as
#   "differentially transcribed". 
#
# ** Downsampling is performed with samtools, which works by hashing
# the template names (thus keeping pair-reads). Note that this may 
# reduce the size of the dataset!
#
# ** The final BAM files are again downsampled between conditions
# to introduce some noise, since they are generated from the same file.
# The count data will reflect this, but not the different "lists",
# however this is easy to adjust afterwards. It appears that templates 
# are affected differently between RFP and RNA by samtools.
######################################################################

__version__ = "0.2"

import os
# set path relative to current directory
thisPath = os.path.dirname(os.path.realpath(__file__))

import shlex
import logging
import subprocess as sp

import pysam
import numpy as np
import pandas as pd


def parse_options():
    
    import argparse as argp
    
    # set parser
    parser = argp.ArgumentParser(description="""Simulates 'differential translational'
        data by downsampling selected subsets of regions from existing BAM files.
        The script generates two conditions from the existing data.""")
    parser.add_argument('-iRFP', '--inputRFP', dest='inputBamRFP', help="""List of input
        RFP BAM file(s), required. Each file will be treated as a replicate from the 
        same condition, in particular subsampling is based on reference names of one of 
        the files, assuming that they all have the same '@SQ SN' header. The order of the 
        files need to match for RFP and RNA if more than one file is given.""", 
        nargs='+', required=True)
    parser.add_argument('-iRNA', '--inputRNA', dest='inputBamRNA', help="""List of input
        RNA BAM file(s), required. Each file will be treated as a replicate from the 
        same condition, in particular subsampling is based on reference names of one of 
        the files, assuming that they all have the same '@SQ SN' header. The order of the 
        files need to match for RFP and RNA if more than one file is given.""", 
        nargs='+', required=True)
    # GROUP 1: TRANSCRIPTIONAL REGULATION, RNA + RFP
    parser.add_argument('--up', help="""Proportion of regions/genes to be 
        'upregulated' overall, for both RNA and RFP. These are concordant
        changes in transcriptional regulation, required.""", type=float, required=True)
    parser.add_argument('--upF', help="""Fraction of templates/pairs to be 
        kept for 'upregulated' genes, required.""", type=float, required=True)
    parser.add_argument('--down', help="""Proportion of regions/genes to be 
        'downregulated' overall, for both RNA and RFP. These are concordant
        changes in transcriptional regulation, required.""", type=float, required=True)
    parser.add_argument('--downF', help="""Fraction of templates/pairs 
        to be kept for 'downregulated' genes, required.""", type=float, required=True)
    # GROUP 2: TRANSLATIONAL REGULATION, RFP ONLY
    parser.add_argument('--rfp-up', help="""Proportion of regions/genes to be 
        ''differentially translated (up)', which are exclusive to RFP, required.""",
        type=float, required=True)
    parser.add_argument('--rfp-upF', help="""Fraction of templates/pairs to be 
        kept for 'upregulated' genes, required.""", type=float, required=True)
    parser.add_argument('--rfp-down', help="""Proportion of regions/genes to be 
        ''differentially translated (down)', which are exclusive to RFP, required.""",
        type=float, required=True)
    parser.add_argument('--rfp-downF', help="""Fraction of templates/pairs 
        to be kept for 'downregulated' genes, required.""", type=float, required=True)
    # GROUP 3: "BUFFERING", RNA ONLY
    parser.add_argument('--rna-up', help="""Proportion of regions/genes to be 
        ''differentially translated (up)', which are exclusive to RNA, required.""",
        type=float, required=True)
    parser.add_argument('--rna-upF', help="""Fraction of templates/pairs to be 
        kept for 'upregulated' genes, required.""", type=float, required=True)
    parser.add_argument('--rna-down', help="""Proportion of regions/genes to be 
        ''differentially translated (down)', which are exclusive to RNA, required.""",
        type=float, required=True)
    parser.add_argument('--rna-downF', help="""Fraction of templates/pairs 
        to be kept for 'downregulated' genes, required.""", type=float, required=True)
    
    parser.add_argument('-d', '--description', dest='sfx', 
                        help='Simulation case description')
    parser.add_argument('--insert', help="""Where to insert the subscript for RFP files
                        which are length-offset adjusted.""", type=str)
    parser.add_argument('--log', help='Log file.')
    parser.add_argument('--version', action='version', version=__version__)
    
    args = parser.parse_args()
     
    # get required arguments if missing
    if not args.sfx:
        args.sfx = 'sim'
    if not args.log:
        logFile = thisPath + '/' + __file__[:-2] + 'log'
        try:
            os.remove(logFile)
        except OSError:
            pass
        args.log = logFile
     
    return args


def extract_reads_from_multiple_regions(inputBam, outputBam, regsList):
    """ Extract all reads that map to a list of regions (regions 
        must be present) and create new bam file.
        
        *** Currently header is simply copied over to the new file, 
        however this includes all other regions that are not
        represented. These files are not kept on their own...
        
        ** Call to samtools handled via subprocess...
        
        Arguments
        ---------
        inputBam: string
            The path/filename of the input bam file
        outputBam: string
            The path/filename of the new bam file.
            The file is overwritten if it exists.
        regsList: string
            The path/filename containing the list of 
            regions to extract from inputBam
    """
    
    if os.path.exists(outputBam):
        os.remove(outputBam)
    filen = os.path.splitext(outputBam)
    outputSam = outputBam + ".sam"
    outputBam += ".bam"
    inputBam += ".bam"
    
    # cannot except CalledProcessError with shell=False...?
    with open(outputSam, 'a') as f:
        cmd = '''samtools view -H ''' + inputBam
        p = sp.call(shlex.split(cmd), stdout=f, shell=False)
        cmd = '''xargs -a ''' + regsList + ''' -I {} samtools view ''' + \
            inputBam + ''' {} '''
        p = sp.call(shlex.split(cmd), stdout=f, shell=False)
    # [W::bgzf_read_block] EOF marker is absent. The input is probably truncated ??
    cmd = '''samtools view -b -h -o ''' + outputBam + """ """ + outputSam
    p = sp.call(shlex.split(cmd), shell=False)
    cmd = ''' rm ''' + outputSam
    p = sp.call(shlex.split(cmd), shell=False)

        
def subsample_bam_file(s, p, inputBam, outputBam):
    """ Subsample BAM files.

        Arguments
        ---------
        s: int
            seed
        p: int
            fraction of templates/pairs to be kept
        inputBam: string
            The path/filename of the input bam file
        outputBam: string
            The path/filename of the new bam file.
    """
    
    s, p = int(s), float(p)
    seed = str(s+p)

    inputBam += ".bam"
    with open(outputBam+".bam", 'w') as f:
        cmd = '''samtools view -b -s ''' + seed + ''' ''' + inputBam
        p = sp.call(shlex.split(cmd), stdout=f, shell=False)
    
    
def merge_bam_file(files2Merge, headerFile, outputBam):
    """ Merge BAM files, sort and index.
    
        *** Does not currently handle different options
        avaiable in samtools/pysam, for instance the
        original header is copied over to the resulting file,
        so this must be given.

        Arguments
        ---------
        files2Merge: list
            List of BAM files to be merged
        headerFile: str
            Original file name if header is to be copied over  
            ** required but no default/no check
        outputBam: string
            The path/filename of the new bam file.
    """
    
    files2Merge = [ str(filen) + ".bam" for filen in files2Merge ]
    headerFile += ".bam"
    outputBamUnsorted = outputBam + ".unsorted.bam"
    outputBam += ".bam"
    
    cmd = '''samtools merge -h ''' + headerFile + ''' ''' + outputBamUnsorted
    for filen in files2Merge:
        cmd += ''' ''' + filen
    p = sp.call(shlex.split(cmd), shell=False)
    cmd = '''samtools sort ''' + outputBamUnsorted + ''' -@10 -o ''' + outputBam
    p = sp.call(shlex.split(cmd), shell=False)
    cmd = '''samtools index ''' + outputBam
    p = sp.call(shlex.split(cmd), shell=False)
    
    
def main():
    """ 
        Input parameters (main -h or parse_options)
        ----------------

        Returns
        -------
            Writes to file region/gene lists corresponding to
            RFP up/RNA down (translational induction) and all others.
            Also writes data (redundantly) in format suitable for visualisation.
            Creates the new BAM files (and a lot of temporary files!)
    """
    
    args = parse_options()
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(args.log)
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s %(asctime)s - %(name)s : %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    msg = """Using following parameters for GROUP 1 (RNA + RFP): 
        up={}, upF={}, down={}, downF={}""".format(args.up,args.upF,
                                                    args.down,args.downF)
    logger.info(msg)
    msg = """Using following parameters for GROUP 2 (RFP ONLY): 
        up={}, upF={}, down={}, downF={}""".format(args.rfp_up,args.rfp_upF,
                                                   args.rfp_down,args.rfp_downF)
    logger.info(msg)
    msg = """Using following parameters for GROUP 3 (RNA ONLY): 
        up={}, upF={}, down={}, downF={}""".format(args.rna_up,args.rna_upF,
                                                   args.rna_down,args.rna_downF)
    logger.info(msg)
    
    missingInputFiles, notBamExt = [], []
    RFPwoExt, RNAwoExt = [], []
    for inputFiles in zip(args.inputBamRFP, args.inputBamRNA):
        for i, files in enumerate(inputFiles):
            if not os.path.exists(files):
                missingInputFiles.append(files)
            filen = os.path.splitext(files)
            if filen[1][1:] != 'bam':
                notBamExt.append(files)
            else:
                if i == 0:
                    RFPwoExt.append(filen[0])
                else:
                    RNAwoExt.append(filen[0])
    
    inputBams = {}
    inputBams['RFP'], inputBams['RNA'] = RFPwoExt, RNAwoExt
    nReplicates = len(RFPwoExt)
                
    if len(missingInputFiles) > 0:
        msg = """Missing input file(s): {}. 
               Terminating.\n""".format(missingInputFiles)
        raise ValueError(msg)

    if len(notBamExt) > 0:
        msg = """Input file(s): {} is/are not BAM file(s)! 
               Terminating.\n""".format(notBamExt)
        raise ValueError(msg)
    
    # Get all reference names and mapped/unmapped read counts.
    logger.info("Get reference names and counts:")
    idxstatsAll = {}
    for dat, inputFiles in inputBams.items():
        idxDf = pd.DataFrame()
        for rep, f in enumerate(inputFiles):
            msg = "{}: {} - {}".format(dat,rep+1, os.path.split(f)[1])
            logger.info(msg)
            # pysam.idxstats returns a long string...!
            idxstats = pysam.idxstats(f + ".bam").splitlines() 
            lines = [ line.strip().split('\t') for line in idxstats]
            idx = [ line[0] for line in lines ]
            cols = [ line[2] for line in lines ]
            colName = "mapped_reads_" + str(rep+1)
            df = pd.DataFrame(cols, index=idx, columns=[colName], dtype=int)
            idxDf = pd.concat([idxDf,df], axis=1)
        idxstatsAll[dat] = idxDf
    
    # remove unmapped reads, we will not be using these if any
    try:
        idxstatsAll['RFP'].drop('*', inplace=True)
        idxstatsAll['RNA'].drop('*', inplace=True)
    except ValueError:
        pass
    
    # Save to file.
    idxstatsAll['RFP'].index.rename('sequence_name', inplace=True)
    dfName = "mapped_read_counts_RFP.csv"
    idxstatsAll['RFP'].to_csv(dfName)
    idxstatsAll['RNA'].index.rename('sequence_name', inplace=True)
    dfName = "mapped_read_counts_RNA.csv"
    idxstatsAll['RNA'].to_csv(dfName)
    
    # Select the regions (genes/reference sequence names).
    # The simplest way to do this is to use one of the replicates, since
    # anyway this is done more or less randomly. 
    # The '@SQ SN' must have the same reference sequence names, which 
    # should be the case if they are from the same experiment.
    rep = "mapped_reads_1"
    
    # We do not want to downsample very low read counts.
    q3RFP = idxstatsAll['RFP'][rep].quantile(0.3)
    q3RNA = idxstatsAll['RNA'][rep].quantile(0.3)
    
    condition = (idxstatsAll['RFP'][rep] > q3RFP)  & \
        (idxstatsAll['RNA'][rep] > q3RNA)
    selection = idxstatsAll['RFP'][rep][condition]
    
    # The number of regions in each group is a proportion
    # of all regions, but they are actually subsampled from
    # this selection.
    allRegions = len(idxstatsAll['RFP'])
    # number of regions group 1
    up = int(args.up*allRegions)
    down = int(args.down*allRegions)
    group1 = up + down
    # number of regions group 2
    rfpUp = int(args.rfp_up*allRegions)
    rfpDown = int(args.rfp_down*allRegions)
    group2 = rfpUp + rfpDown
    # number of regions group 3
    rnaUp = int(args.rna_up*allRegions)
    rnaDown = int(args.rna_down*allRegions)
    group3 = rnaUp + rnaDown
    
    msg = """Number of regions in each group: 
        GROUP 1={}, GROUP 2={}, GROUP 3={}""".format(group1,group2,group3)
    logger.info(msg)
    
    allGroups = group1 + group2 + group3
    if allGroups > len(selection):
        msg = """Cannot sample a total of {} regions from {}! 
               Change threshold(s) for subsetting 
               data.\n""".format(allGroups,len(selection))
        raise ValueError(msg)
    
    # Subset regions to form groups, and remove them 
    # from selection. There is no overlap between groups.
    # Keep track of all the selected regions.
    # For each condition, we also create lists of unchanged regions
    # to facilitate extraction.
    seed = np.random.randint(100, size=1)[0]
    upRegions = selection.sample(n=up, random_state=seed)
    selection = selection.drop(selection.loc[upRegions.index].index)
    selectedRegions = upRegions.index
    cond2RegionsRFP = upRegions.index
    cond2RegionsRNA = upRegions.index
    
    # Repeat for regions that will be "down-regulated" ...
    seed = np.random.randint(100, size=1)[0]
    downRegions = selection.sample(n=down, random_state=seed)
    selection = selection.drop(selection.loc[downRegions.index].index)
    selectedRegions = selectedRegions.append(downRegions.index)
    cond1RegionsRFP = downRegions.index
    cond1RegionsRNA = downRegions.index
    
    # ... "differentially translated" (group 2)
    seed = np.random.randint(100, size=1)[0]
    rfpUpRegions = selection.sample(n=rfpUp, random_state=seed)
    selection = selection.drop(selection.loc[rfpUpRegions.index].index)
    selectedRegions = selectedRegions.append(rfpUpRegions.index)
    seed = np.random.randint(100, size=1)[0]
    rfpDownRegions = selection.sample(n=rfpDown, random_state=seed)
    selection = selection.drop(selection.loc[rfpDownRegions.index].index)
    selectedRegions = selectedRegions.append(rfpDownRegions.index)
    
    cond1RegionsRFP = cond1RegionsRFP.append(rfpDownRegions.index)
    cond2RegionsRFP = cond2RegionsRFP.append(rfpUpRegions.index)
    
    cond1RegionsRNA = cond1RegionsRNA.append(rfpUpRegions.index)
    cond2RegionsRNA = cond2RegionsRNA.append(rfpUpRegions.index)
    cond1RegionsRNA = cond1RegionsRNA.append(rfpDownRegions.index)
    cond2RegionsRNA = cond2RegionsRNA.append(rfpDownRegions.index)
    
    # ... "buffered" (group 3)
    seed = np.random.randint(100, size=1)[0]
    rnaUpRegions = selection.sample(n=rnaUp, random_state=seed)
    selection = selection.drop(selection.loc[rnaUpRegions.index].index)
    selectedRegions = selectedRegions.append(rnaUpRegions.index)
    seed = np.random.randint(100, size=1)[0]
    rnaDownRegions = selection.sample(n=rnaDown, random_state=seed)
    selection = selection.drop(selection.loc[rnaDownRegions.index].index)
    selectedRegions = selectedRegions.append(rnaDownRegions.index)
    
    cond1RegionsRNA = cond1RegionsRNA.append(rnaDownRegions.index)
    cond2RegionsRNA = cond2RegionsRNA.append(rnaUpRegions.index)
    
    cond1RegionsRFP = cond1RegionsRFP.append(rnaUpRegions.index)
    cond2RegionsRFP = cond2RegionsRFP.append(rnaUpRegions.index)
    cond1RegionsRFP = cond1RegionsRFP.append(rnaDownRegions.index)
    cond2RegionsRFP = cond2RegionsRFP.append(rnaDownRegions.index)
    
    # Get all regions that will remain unchanged.
    unchangedRegions = idxstatsAll['RFP'].drop(selectedRegions).index
    cond1RegionsRFP = cond1RegionsRFP.append(unchangedRegions)
    cond2RegionsRFP = cond2RegionsRFP.append(unchangedRegions)
    cond1RegionsRNA = cond1RegionsRNA.append(unchangedRegions)
    cond2RegionsRNA = cond2RegionsRNA.append(unchangedRegions)

    msg = """Total number of unchanged regions: {}""".format(len(unchangedRegions))
    logger.info(msg)
    
    # Only checking to make sure...
    if selectedRegions.get_duplicates() or \
        unchangedRegions.get_duplicates():
        msg = """Duplicate indices...!\n"""
        raise ValueError(msg)
    if not unchangedRegions.intersection(selectedRegions).empty:
        msg = """Overlap between selected and unchanged regions...!\n"""
        raise ValueError(msg)

    # Keep track... we need these to tell samtools which regions to 
    # extract from our original files. 
    listTracking = {'up':'rfp-rna-up', 'down':'rfp-rna-down', 
                    'rfp-up':'rfp-up', 'rfp-down':'rfp-down', 'rna-up':'rna-up',
                    'rna-down':'rna-down', 'unchanged':'unchanged',
                    'cond1rfp':'rfp-unchanged-cond1', 'cond2rfp':'rfp-unchanged-cond2',
                    'cond1rna':'rna-unchanged-cond1', 'cond2rna':'rna-unchanged-cond2',}
    np.savetxt(listTracking['up'], upRegions.index, delimiter=",", fmt='%s')
    np.savetxt(listTracking['down'], downRegions.index, delimiter=",", fmt='%s')
    np.savetxt(listTracking['rfp-up'], rfpUpRegions.index, delimiter=",", fmt='%s')
    np.savetxt(listTracking['rfp-down'], rfpDownRegions.index, delimiter=",", fmt='%s')
    np.savetxt(listTracking['rna-up'], rnaUpRegions.index, delimiter=",", fmt='%s')
    np.savetxt(listTracking['rna-down'], rnaDownRegions.index, delimiter=",", fmt='%s')
    np.savetxt(listTracking['unchanged'], unchangedRegions, delimiter=",", fmt='%s')
    np.savetxt(listTracking['cond1rfp'], cond1RegionsRFP, delimiter=",", fmt='%s')
    np.savetxt(listTracking['cond2rfp'], cond2RegionsRFP, delimiter=",", fmt='%s')
    np.savetxt(listTracking['cond1rna'], cond1RegionsRNA, delimiter=",", fmt='%s')
    np.savetxt(listTracking['cond2rna'], cond2RegionsRNA, delimiter=",", fmt='%s')
    
    # Sort and generate all BAM files (and intermediate SAM files!)
    msg = """Sort and generate all BAM files..."""
    logger.info(msg)
    seeds = np.random.randint(100, size=2*nReplicates).reshape((2,nReplicates))
    
    fileTracking = { key:{'cond1':[],'cond2':[]} for key in ['RFP', 'RNA']}
    cond1 = fileTracking['RFP']['cond1']
    cond2 = fileTracking['RFP']['cond2']
    outputBams = { key:{'cond1':[],'cond2':[]} for key in ['RFP', 'RNA']}
    out1 = outputBams['RFP']['cond1']
    out2 = outputBams['RFP']['cond2']
    for rep, inputFile in enumerate(inputBams['RFP']):
        msg = """Using input: {}""".format(inputFile)
        logger.info(msg)
        # condition 1
        suffix = "." + args.sfx + ".cond1"
        if args.insert:
            f = inputFile.split(args.insert)
            outFile = f[0] + suffix[1:] + '.' + args.insert + f[1]
        else:
            outFile = inputFile + suffix
        out1.append(outFile)
        f1 = outFile + "." + listTracking['up'] + ".sub"
        f2 = outFile + "." + listTracking['rfp-up'] + ".sub"
        f3 = outFile + "." + listTracking['cond1rfp']
        cond1.append((f1, f2, f3))
        extract_reads_from_multiple_regions(inputFile, f1[:-4], listTracking['up'])
        extract_reads_from_multiple_regions(inputFile, f2[:-4], listTracking['rfp-up'])
        extract_reads_from_multiple_regions(inputFile, f3, listTracking['cond1rfp'])
        subsample_bam_file(seeds[0,rep], args.upF, f1[:-4], f1)
        subsample_bam_file(np.random.randint(100, size=1)[0], args.rfp_upF, f2[:-4], f2)
        
        # condition 2
        suffix = "." + args.sfx + ".cond2"
        if args.insert:
            f = inputFile.split(args.insert)
            outFile = f[0] + suffix[1:] + '.' + args.insert + f[1]
        else:
            outFile = inputFile + suffix
        out2.append(outFile)
        f1 = outFile + "." + listTracking['down'] + ".sub"
        f2 = outFile + "." + listTracking['rfp-down'] + ".sub"
        f3 = outFile + "." + listTracking['cond2rfp']
        cond2.append((f1, f2, f3))
        extract_reads_from_multiple_regions(inputFile, f1[:-4], listTracking['down'])
        extract_reads_from_multiple_regions(inputFile, f2[:-4], listTracking['rfp-down'])
        extract_reads_from_multiple_regions(inputFile, f3, listTracking['cond2rfp'])
        subsample_bam_file(seeds[1,rep], args.downF, f1[:-4], f1)
        subsample_bam_file(np.random.randint(100, size=1)[0], args.rfp_downF, f2[:-4], f2)
        
        msg = """These files will be created: 
            cond1:{}
            cond2:{}""".format(out1[-1],out2[-1])
        logger.info(msg)
        
    cond1 = fileTracking['RNA']['cond1']
    cond2 = fileTracking['RNA']['cond2']
    out1 = outputBams['RNA']['cond1']
    out2 = outputBams['RNA']['cond2']
    for rep, inputFile in enumerate(inputBams['RNA']):
        msg = """Using input: {}""".format(inputFile)
        logger.info(msg)
        # condition 1
        suffix = "." + args.sfx + ".cond1"
        outFile = inputFile + suffix
        out1.append(outFile)
        f1 = outFile + "." + listTracking['up'] + ".sub"
        f2 = outFile + "." + listTracking['rna-up'] + ".sub"
        f3 = outFile + "." + listTracking['cond1rna']
        cond1.append((f1, f2, f3))
        extract_reads_from_multiple_regions(inputFile, f1[:-4], listTracking['up'])
        extract_reads_from_multiple_regions(inputFile, f2[:-4], listTracking['rna-up'])
        extract_reads_from_multiple_regions(inputFile, f3, listTracking['cond1rna'])
        subsample_bam_file(seeds[0,rep], args.upF, f1[:-4], f1)
        subsample_bam_file(np.random.randint(100, size=1)[0], args.rna_upF, f2[:-4], f2)
        
        # condition 2
        suffix = "." + args.sfx + ".cond2"
        outFile = inputFile + suffix
        out2.append(outFile)
        f1 = outFile + "." + listTracking['down'] + ".sub"
        f2 = outFile + "." + listTracking['rna-down'] + ".sub"
        f3 = outFile + "." + listTracking['cond2rna']
        cond2.append((f1, f2, f3))
        extract_reads_from_multiple_regions(inputFile, f1[:-4], listTracking['down'])
        extract_reads_from_multiple_regions(inputFile, f2[:-4], listTracking['rna-down'])
        extract_reads_from_multiple_regions(inputFile, f3, listTracking['cond2rna'])
        subsample_bam_file(seeds[1,rep], args.downF, f1[:-4], f1)
        subsample_bam_file(np.random.randint(100, size=1)[0], args.rna_downF, f2[:-4], f2)
        
        msg = """These files will be created: 
            cond1:{}
            cond2:{}""".format(out1[-1],out2[-1])
        logger.info(msg)
        
    # Merge all reads with original headers ** header should remain unchanged...
    # Sort and index
    msg = """Merging all files..."""
    logger.info(msg)
    for dat, conditions in fileTracking.items():
        for cond, files in conditions.items():
            for rep, repFiles in enumerate(files):
                files2Merge = list(repFiles)
                inputBam = inputBams[dat][rep]
                outputBam = outputBams[dat][cond][rep]
                msg="""{}, {}, using {} to generate {}""".format(dat,cond,repFiles,outputBam)
                logger.info(msg)
                merge_bam_file(files2Merge, inputBam, outputBam)
    
    # Introduce some noise in the data.
    seeds1 = {'cond1':[],'cond2':[]}
    seeds2 = {'cond1':[],'cond2':[]} 
    for cond in ['cond1', 'cond2']:
        for rep in range(nReplicates):
           seeds1[cond].append(np.random.randint(100, size=1)[0])
           seeds2[cond].append(np.random.randint(100, size=1)[0])

    for dat, conditions in outputBams.items():
        for cond, outputFiles in conditions.items():
            for rep, repFiles in enumerate(outputFiles):
                tmpf = repFiles + ".bak.bam"
                os.rename(repFiles + ".bam", tmpf)
                subsample_bam_file(seeds1[cond][rep], .87, tmpf[:-4], tmpf[:-4] + ".bak")
                subsample_bam_file(seeds2[cond][rep], .82, tmpf[:-4] + ".bak", repFiles)
                os.remove(tmpf) 
                os.remove(tmpf[:-4] + ".bak.bam")
                tmpf = repFiles + '.bam'
                cmd = '''samtools index ''' + tmpf
                p = sp.call(shlex.split(cmd), shell=False)

    # Get counts for new data.
    idxstatsNew = {}
    for dat, conditions in outputBams.items():
        idxDf = pd.DataFrame()
        for cond, outputFiles in conditions.items():
            for rep, repFiles in enumerate(outputFiles):
                msg = "New data {}, {}: {} - {}".format(dat,cond,rep+1,os.path.split(repFiles)[1])
                logger.info(msg)
                # pysam.idxstats returns a long string...!
                idxstats = pysam.idxstats(repFiles + ".bam").splitlines() 
                lines = [ line.strip().split('\t') for line in idxstats]
                idx = [ line[0] for line in lines ]
                cols = [ line[2] for line in lines ]
                colName = "mapped_reads_" + str(cond) + '_' + str(rep+1)
                df = pd.DataFrame(cols, index=idx, columns=[colName], dtype=int)
                idxDf = pd.concat([idxDf,df], axis=1)
        idxstatsNew[dat] = idxDf
        
    # Save to file.
    idxstatsNew['RFP'].index.rename('sequence_name', inplace=True)
    dfName = "mapped_read_counts_RFP_" + args.sfx + ".csv"
    idxstatsNew['RFP'].to_csv(dfName)
    idxstatsNew['RNA'].index.rename('sequence_name', inplace=True)
    dfName = "mapped_read_counts_RNA_" + args.sfx + ".csv"
    idxstatsNew['RNA'].to_csv(dfName)
    
    # remove files
    for dat, conditions in fileTracking.items():
        for cond, files in conditions.items():
            for rep, repFiles in enumerate(files):
                for bygone in repFiles:
                    try:
                        filen = bygone + '.bam'
                        os.remove(filen)
                        msg = "Removing: {}".format(filen)
                        logger.info(msg)
                        filen = os.path.splitext(bygone)
                        if filen[1][1:] == 'sub':
                            filen = filen[0] + '.bam'
                            os.remove(filen)
                            msg = "Removing: {}".format(filen)
                            logger.info(msg)
                    except OSError:
                        msg = "Could not remove: {}".format(bygone)
                        logger.info(msg)
    import glob
    for filen in glob.glob(thisPath+'/*unsorted*'):
        try:
            os.remove(filen) 
            msg = "Removing: {}".format(filen)
            logger.info(msg)
        except OSError:
            msg = "Could not remove: {}".format(filen)
            logger.info(msg)


if __name__ == "__main__":
     
    main()
