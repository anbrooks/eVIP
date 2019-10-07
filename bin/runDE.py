########################################################################
# Adapted from runDE.py in FLAIR
# https://github.com/BrooksLabUCSC/flair
# Original author: Cameron M. Soulette
#
########################################################################

import os, sys
import pandas as pd
import numpy as np

from rpy2 import robjects
from rpy2.robjects import r,pandas2ri, Formula
from rpy2.robjects.lib import grid
pandas2ri.activate()
R = robjects.r

#supressing rpy2 warnings
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)


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
def main():
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
    formulaDF     = pd.read_csv(formula,header=0, sep="\t",index_col=0)
    print("pandas df")
    print(formulaDF)

    samples =  formulaDF.index.tolist()
    R.assign('samples',samples)

    sampleTable = pandas2ri.py2ri(formulaDF)
    R.assign('sampleTable',sampleTable)

    R('colnames(sampleTable) = c("condition")')

    R('sampleTable$ID <- seq.int(nrow(sampleTable))')


    # R('sampleTable <- data.frame(samples=samples, condition=condition)')


    print("R df:")
    R('print(sampleTable)')


    #locate kallisto files
    R('files <- file.path(indir, samples, "abundance.h5")')
    R('all(file.exists(files))')

    #tximport conversion to gene
    R('txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)')
    R('rownames(sampleTable) <- samples')
    R('print(sampleTable)')

    #DESeq
    R('dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)')


    R('colData(dds)$condition<-factor(colData(dds)$condition, levels=c("group1","group2"))')

    print("~~~~~colData")
    R('print(colData(dds))')

    R('GFP_659_dds<-DESeq(dds)')
    R('GFP_659_res<-results(GFP_659_dds)')
    R('GFP_659_res<-GFP_659_res[order(GFP_659_res$padj),]')

    R('print(head(GFP_659_res))')

    #writing deseq2 results to a file
    # R(write.csv(as.data.frame(GFP_659_res),file='GFP_vs_659_results_deseq2.csv'))


# ## Merge with normalized count data
# GFP_659_resdata <- merge(as.data.frame(GFP_659_res), as.data.frame(counts(GFP_659_dds, normalized=TRUE)), by="row.names", sort=FALSE)
# names(GFP_659_resdata)[1] <- "Gene"
# head(GFP_659_resdata)
# mcols(GFP_659_res,use.names=TRUE)
# write.csv(as.data.frame(GFP_659_res),file='GFP_vs_659_results_deseq2.csv')
# write.csv(as.data.frame(GFP_659_resdata),file='GFP_vs_659_results_deseq2_counts.csv')



    """

    # make the quant DF
    quantDF  = pd.read_table(matrix, header=0, sep='\t', index_col=0)
    df = pandas2ri.py2ri(quantDF)



    ### RUN DESEQ2 ###
    R.assign('df', df)
    R.assign('sampleTable', sampleTable)
    R.assign('design',design)
    R('dds <- DESeqDataSetFromMatrix(countData = df, colData = sampleTable, design = design)')
    R('dds <- DESeq(dds)')
    R('name <- grep("condition", resultsNames(dds), value=TRUE)')

    ###
    ###
    # Get Results and shrinkage values
    res    = R('results(dds, name=name)')
    resLFC = R('lfcShrink(dds, coef=name)')
    vsd    = R('vst(dds,blind=FALSE)')
    resdf  = robjects.r['as.data.frame'](res)
    reslfc = robjects.r['as.data.frame'](resLFC)
    dds    = R('dds')



    data_folder = os.path.join(os.getcwd(), outdir)
    lfcOut = os.path.join(data_folder, "%s_%s_v_%s_deseq2_results_shrinkage.tsv"  % (prefix,group1,group2))
    resOut = os.path.join(data_folder, "%s_%s_v_%s_deseq2_results.tsv"  % (prefix,group1,group2))

    robjects.r['write.table'](reslfc, file=lfcOut, quote=False, sep="\t")
    robjects.r['write.table'](resdf, file=resOut, quote=False, sep="\t")

    """
if __name__ == "__main__":
    main()
