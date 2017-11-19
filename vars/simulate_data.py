######################################################################
# This file is an ad hoc script used to simulate data for the b-tea
# pipeline. TE can be estimated using a simulation approach based on 
# modifying real data to introduce differences in RNA/RFP data. 
# This script uses samtools to subsample existing RNA/RFP BAM files, 
# except for a proportion of known regions/genes, that "should" show 
# up as differentially translated. 
#
# ** Downsampling is performed with samtools, which works by hashing
# the template names (thus keeping pair-reads). Note that this may 
# significantly reduce the size of the dataset!
#
# ** The script currently generates many temporary files: these are not
# deleted automatically!  
######################################################################

__version__ = "0.1"

import os
# set path relative to current directory
thisPath = os.path.dirname(os.path.realpath(__file__))

import shlex
import subprocess as sp

import pysam
import numpy as np
import pandas as pd


def parse_options(activeOptions):
    """ Parse options/arguments.
        No user input if active options are not set.
        
        Arguments:
        ----------
        activeOptions: list
            The list of active options
        
        Returns:
        --------
        args: dict
            The active options
    """
    import argparse as argp

    # set parser
    parser = argp.ArgumentParser()
     
    for opts in activeOptions:
        if opts == 'iRFP':
            parser.add_argument('-iRFP', '--inputRFP', dest='inputBamRFP',
                                help='Input RFP BAM file, required',
                                required=True)
        elif opts == 'iRNA':
            parser.add_argument('-iRNA', '--inputRNA', dest='inputBamRNA', 
                                help='Input RNA BAM file, required', 
                                required=True)
        elif opts == 'pcTrl':
            msg = """ Percentage of regions/genes untouched (RFP), required. 
                      These are downsampled for RNA only. """
            parser.add_argument('-pcTrl', '--pcTrl', dest='pcTrl', 
                                help=msg, type=int, required=True)
        elif opts == 'fracAll':
            msg = """ Percentage of overall templates/pairs to throw away, required """
            parser.add_argument('-fracAll', '--fracAll', dest='pcFracAll', 
                                help=msg, type=int, required=True)
        elif opts == 'fracRNA':
            msg = """ Percentage of RNA templates/pairs to throw away, required """
            parser.add_argument('-fracRNA', '--fracRNA', dest='pcFracRNA', 
                                help=msg, type=int, required=True)
        elif opts == 'nR':
            parser.add_argument('-nR', '--nReplicates', dest='nReplicates', 
                                help='Number of replicates to genmerate.', type=int, 
                                default = 3)
        elif opts == 'd':
            parser.add_argument('-d', '--description', dest='sfx', 
                                help='Simulation case description')
        else:
            print('Warning: option {0} not set, ignoring... '.format(opts))
         
    parser.add_argument('--version', action='version', version=__version__)
    args = vars(parser.parse_args()) # dictionary
     
    # get required arguments if missing
    if not args.get('sfx'):
        args['sfx'] = 'sim'
     
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
    
    s, p = int(s), int(p)
    if p < 10: p = "0" + str(p)
    seed = str(str(s) + "." + str(p))

    inputBam += ".bam"
    with open(outputBam+".bam", 'w') as f:
        cmd = '''samtools view -b -s ''' + seed + ''' ''' + inputBam
        p = sp.call(shlex.split(cmd), stdout=f, shell=False)
    
    
def merge_bam_file(files2Merge, headerFile, outputBam):
    """ Merge BAM files, sort and index.
    
        *** Does not currently handle different options
        avaiable in sammtools/pysam, for instance the
        original header is copied over to the resulting file,
        so this must be given.

        Arguments
        ---------
        files2Merge: list
            List of BAM files to be merged
        headerFile: str
            Original file name if header is to be copied over  
            ** required but no defaukt/no check
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
    cmd = '''samtools sort ''' + outputBamUnsorted + ''' -o ''' + outputBam
    p = sp.call(shlex.split(cmd), shell=False)
    cmd = '''samtools index ''' + outputBam
    p = sp.call(shlex.split(cmd), shell=False)
    
    
def main():
    """ 
        Input parameters (main -h or parse_options)
        ----------------
            iRFP
            iRNA
            pcTrl
            fracAll
            fracRNA
            nR
            d
        
        Returns
        -------
            Writes to file region/gene lists corresponding to
            RFP up/RNA down (translational induction) and all others.
            Also writes data (redundantly) in format suitable for visualisation.
            Creates the new BAM files (and a lot of temporary files!)
    """
    
    activeOptions = parse_options(['iRFP','iRNA','pcTrl','fracAll','fracRNA','nR','d'])
    
    inputBamRFP = activeOptions['inputBamRFP']
    inputBamRNA = activeOptions['inputBamRNA']
    pcTrl = activeOptions['pcTrl']
    pcFracAll = activeOptions['pcFracAll']
    pcFracRNA = activeOptions['pcFracRNA']
    nReplicates = activeOptions['nReplicates']
    sfx = activeOptions['sfx']
    
    # Set proportion (percentage) of templates/pairs to be kept for downsampling.
    pcFracRNA = 100 - pcFracRNA # this is where we want to see "translational induction" 
    pcFracAll = 100 - pcFracAll
    
    missingInputFiles, notBamExt, inputFilesWoExt = [], [], []
    for inputFiles in [inputBamRFP, inputBamRNA]:
        if not os.path.exists(inputFiles):
            missingInputFiles.append(missingInputFiles)
        filen = os.path.splitext(inputFiles)
        if filen[1][1:] != 'bam':
            notBamExt.append(inputFiles)
        else:
            inputFilesWoExt.append(filen[0])
            
    if len(missingInputFiles) > 0:
        msg = """Missing input file(s): {}. 
               Terminating.\n""".format(missingInputFiles)
        raise ValueError(msg)

    if len(notBamExt) > 0:
        msg = """Some input files: {} are not BAM files! 
               Terminating.\n""".format(notBamExt)
        raise ValueError(msg)

    inputBamRFP, inputBamRNA = inputFilesWoExt[0], inputFilesWoExt[1]

    # Get reference names and counts.
    # pysam.idxstats returns a long string...!
    idxstats = pysam.idxstats(inputBamRNA+".bam").splitlines()
    lines = [ line.strip().split('\t') for line in idxstats]
    idx = [ line[0] for line in lines ]
    cols = [ line[2] for line in lines ]
    idxstatsRNA = pd.DataFrame(cols, index=idx, columns=["mapped_reads"], dtype=int)

    idxstats = pysam.idxstats(inputBamRFP+".bam").splitlines()
    lines = [ line.strip().split('\t') for line in idxstats]
    idx = [ line[0] for line in lines ]
    cols = [ line[2] for line in lines ]
    idxstatsRFP = pd.DataFrame(cols, index=idx, columns=["mapped_reads"], dtype=int)

    nAllRegions = len(idxstatsRFP)
    nSubset = int(pcTrl*nAllRegions*0.01)
    
    # Subsampling is performed above the upper quartile of both RFP/RNA 
    # and excludes upper outliers for RNA
    q3RFP = idxstatsRFP['mapped_reads'].quantile(0.75)
    q3RNA = idxstatsRNA['mapped_reads'].quantile(0.75)
    iqrRNA = q3RNA - idxstatsRNA['mapped_reads'].quantile(0.25)
    fenceHighRNA = q3RNA + 3.*iqrRNA
    
    condition = (idxstatsRFP['mapped_reads'] > q3RFP)  & \
        (idxstatsRNA['mapped_reads'] > q3RNA) & \
            (idxstatsRNA['mapped_reads'] < fenceHighRNA)
    selection = idxstatsRFP[condition]
    if nSubset > len(selection):
        msg = """Cannot sample {} regions from {}! 
               Change threshold for subsetting data.\n""".format(nSubset,len(selection))
        raise ValueError(msg)
    seed = np.random.randint(100, size=1)[0]
    subset = selection.sample(n=nSubset, random_state=seed)
    allOthers = idxstatsRFP.drop(idxstatsRFP.loc[subset.index].index)
    
    # Keep track of files... 
    listTracking = ['trl', 'ctrl'] 
    subset.to_csv(listTracking[0], columns=[], header=False)
    allOthers.to_csv(listTracking[1], columns=[], header=False)
    
    # Sort and generate all BAM files (and intermediate SAM files!)
    replicates = [ "rep" + str(i+1) for i in range(nReplicates) ]
    seedsAll = np.random.randint(100, size=nReplicates)
    seedsRNA = np.random.randint(100, size=nReplicates)
    outputBamRFP, outputBamRNA = [], []
    fileTrackingRFP, fileTrackingRNA = [], []
    
    for i, replicate in enumerate(replicates):
        suffix = "_" + sfx + "." + replicate
        nameRFP = inputBamRFP + suffix
        nameRNA = inputBamRNA + suffix
        outputBamRFP.append(nameRFP)
        outputBamRNA.append(nameRNA)
        
        # RFP
        f1 = nameRFP + "." + listTracking[0]
        f2 = nameRFP + "." + listTracking[1] + ".sub"
        fileTrackingRFP.append((f1, f2))
        extract_reads_from_multiple_regions(inputBamRFP, f1, listTracking[0])
        extract_reads_from_multiple_regions(inputBamRFP, f2[:-4], listTracking[1])
        subsample_bam_file(seedsAll[i], pcFracAll, f2[:-4], f2)
        
        # RNA
        f1 = nameRNA + "." + listTracking[0] + ".sub"
        f2 = nameRNA + "." + listTracking[1] + ".sub"
        fileTrackingRNA.append((f1, f2))
        extract_reads_from_multiple_regions(inputBamRNA, f1[:-4], listTracking[0])
        subsample_bam_file(seedsRNA[i], pcFracRNA, f1[:-4], f1)
        extract_reads_from_multiple_regions(inputBamRNA, f2[:-4], listTracking[1])
        subsample_bam_file(seedsAll[i], pcFracAll, f2[:-4], f2)
            
        
    # Merge all reads with original headers ** header should remain unchanged...
    # Sort and index
    for i, replicate in enumerate(replicates):
        files2Merge = list(fileTrackingRFP[i])
        merge_bam_file(files2Merge, inputBamRFP, outputBamRFP[i])
        files2Merge = list(fileTrackingRNA[i])
        merge_bam_file(files2Merge, inputBamRNA, outputBamRNA[i])
    
    # Get counts for new data.
    idxstatsNewRFP, idxstatsNewRNA = [], []
    for i in range(nReplicates):
        # RFP
        idxstats = pysam.idxstats(outputBamRFP[i]+".bam").splitlines()
        lines = [ line.strip().split('\t') for line in idxstats]
        idx = [ line[0] for line in lines ]
        cols = [ line[2] for line in lines ]
        df = pd.DataFrame(cols, index=idx, columns=["rfp_new"], dtype=int)
        df.sort_index(inplace=True)
        idxstatsNewRFP.append(df)
        # RNA
        idxstats = pysam.idxstats(outputBamRNA[i]+".bam").splitlines()
        lines = [ line.strip().split('\t') for line in idxstats]
        idx = [ line[0] for line in lines ]
        cols = [ line[2] for line in lines ]
        df = pd.DataFrame(cols, index=idx, columns=["rna_new"], dtype=int)
        df.sort_index(inplace=True)
        idxstatsNewRNA.append(df)
        
    # Compare results.
    idxstatsRFP.rename(columns={"mapped_reads":"rfp_old"}, inplace=True)
    idxstatsRNA.rename(columns={"mapped_reads":"rna_old"}, inplace=True)
    for i, replicate in enumerate(replicates):
        dfList = [idxstatsRFP, idxstatsNewRFP[i], idxstatsRNA, idxstatsNewRNA[i]]
        df = pd.concat(dfList, axis=1)
        df['change_rfp'] = df['rfp_old']-df['rfp_new']
        df['change_rna'] = df['rna_old']-df['rna_new']
        
        # 1:  "translational induction" (no changes in RFP counts, but changes in RNA counts)
        # 0:  no changes in both (untouched) incl. zero counts
        # -1: default to downsampled regions RFP/RNA (and RFP/- or -/RNA?)
        # 
        conditions = [ (df['change_rfp'] == 0) & (df['change_rna'] != 0) & \
            (df['rfp_old'] > q3RFP), (df['change_rfp'] == 0) & (df['change_rna'] == 0) ]
        choices = [1, 0]
        df['te_flag'] = np.select(conditions, choices, default=-1)
        trl = df[df['te_flag'] == 1]
        dfName = "mapped_read_counts." + sfx + "." + replicate + ".translated"
        trl.to_csv(dfName, index_label='sequence_name', columns=["rfp_old", "rfp_new", "rna_old", "rna_new"])

        # Now reshape for plotting...
        allRFP = pd.Series(df[['rfp_old','rfp_new']].values.ravel('F'))
        allRNA = pd.Series(df[['rna_old','rna_new']].values.ravel('F'))
        flag = pd.Series(['Original']*len(df)+['Simulated']*len(df))
        
        trl = trl.loc[:,['rfp_new','rna_new']]
        trl.rename(columns={"rfp_new":"RFP", "rna_new":"RNA"}, inplace=True)
        trl.reset_index(drop=True, inplace=True)
        allRFP = allRFP.append(trl['RFP'])
        allRNA = allRNA.append(trl['RNA'])
        flag = flag.append(pd.Series(['Translated']*len(trl)))
        
        allData = pd.DataFrame({'RFP': allRFP, 'RNA':allRNA, 'Data':flag})
        dfName = "scatterAll." + sfx + "." + replicate
        allData.to_csv(dfName, index=False)
        

if __name__ == "__main__":
     
    main()