########################################################################
# Adapted from runDE.py in FLAIR
# https://github.com/BrooksLabUCSC/flair
# Original author: Cameron M. Soulette
#
########################################################################

import os, sys
import pandas as pd
import numpy as np

#supressing rpy2 warnings
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

from rpy2 import robjects
from rpy2.robjects import r,pandas2ri, Formula
from rpy2.robjects.lib import grid
pandas2ri.activate()
R = robjects.r




########################################################################
# CommandLine
########################################################################

class CommandLine(object) :
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    and a standard usage and help,
    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']

    methods:

    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = ' runDE.py - a rpy2 convenience tool to run DESeq2.',
                                             epilog = 'Please feel free to forward any questions/concerns to /dev/null',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s ')
        # Add args
        self.parser.add_argument("--group1"    , action = 'store', required=True,
                                    help='Sample group 1.')
        self.parser.add_argument("--group2"    , action = 'store', required=True,
                                    help='Sample group 2.')
        self.parser.add_argument("--outDir"    , action = 'store', required=True,
                                    help='Write to specified output directory.')
        self.parser.add_argument("--inDir"    , action = 'store', required=True,
                                    help='Directory where kallisto output directories are located')
        self.parser.add_argument("--formula"    , action = 'store', required=True,
                                    help='Formula design matrix.')

        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

# main
def main(group1=None, group2=None, outDir=None, inDir=None, formula=None):
    '''
    main
    '''

    # Command Line Stuff...
    myCommandLine = CommandLine()

    indir      = myCommandLine.args['inDir']
    outdir     = myCommandLine.args['outDir']
    group1     = myCommandLine.args['group1']
    group2     = myCommandLine.args['group2']
    formula    = myCommandLine.args['formula']

    R.assign('indir',indir)
    R.assign('outdir',outdir)

    R.assign('group1',group1)
    R.assign('group2',group2)

    # import
    from rpy2.robjects.packages import importr
    #kallisto processing libraries
    tximportData   = importr('tximportData')
    tximport   = importr('tximport')
    ensembldb   = importr('ensembldb')
    EnsDb_Hsapiens_v86   = importr('EnsDb.Hsapiens.v86')
    #deseq
    methods   = importr('methods')
    deseq     = importr('DESeq2')

    #transcripts to gene, used in tximport
    R('edb <- EnsDb.Hsapiens.v86')
    R('tx2gene = transcripts(edb , columns=c("tx_id", "gene_name"),return.type="DataFrame")')

    # import formula
    formulaDF     = pd.read_csv(formula,header=0, sep="\t")

    samples =  formulaDF.samples.tolist()
    R.assign('samples',samples)

    sampleTable = pandas2ri.py2ri(formulaDF)
    R.assign('sampleTable',sampleTable)

    #locate kallisto files
    R('files <- file.path(indir, samples, "abundance.h5")')
    R('all(file.exists(files))')

    #tximport conversion to gene
    R('txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)')
    R('rownames(sampleTable) <- samples')

    #DESeq
    R('dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)')
    R('colData(dds)$condition<-factor(colData(dds)$condition, levels=c(group1,group2))')

    #printing samples and condition
    R('print(colData(dds))')

    R('dds_<-DESeq(dds)')
    R('res<-results(dds_)')
    R('res<-res[order(res$padj),]')

    print("Deseq2 results")
    R('print(head(res))')

    # writing deseq2 results to a file
    # R(write.csv(as.data.frame(res),file='GFP_vs_659_results_deseq2.csv'))

    Out = os.path.join(outdir, "%s_v_%s_deseq2_results.tsv"  % (group1,group2))
    R.assign('Out',Out)

    R('write.csv(as.data.frame(res),file=Out)')

if __name__ == "__main__":
    main()
