# <Script name>
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb
import math

#############
# CONSTANTS #
#############
DEF_MIN_FOLD_FPKM = 0.0
DEF_MIN_FPKM = 1.0

START_VAL_IDX = 1
INFINITY = 10
#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)


###############
# END CLASSES #
###############
 
########
# MAIN #	
########
def main():
	
    opt_parser = OptionParser()
   
    # Add Options. Required options should have default=None
    opt_parser.add_option("--in_table",
                          dest="in_table",
                          type="string",
                          help="Input table of gene expression values",
                          default=None)
    opt_parser.add_option("--out_table",
                          dest="out_table",
                          type="string",
                          help="""Filtered output table. Filter is based on FPKM
                                  minimum threshold and minimum FPKM difference
                                  cutoffs""",
                          default=None)
    opt_parser.add_option("-i",
                          dest="start_idx",
                          type="int",
                          help="0-based index of the start of the values. Def=%d" % START_VAL_IDX,
                          default=START_VAL_IDX)
    opt_parser.add_option("-l",
                          dest="log2_transform",
                          action="store_true",
                          help="Will log2 transform the data.",
                          default=None)
    opt_parser.add_option("--reformat_gene",
                          dest="gene_col",
                          type="int",
                          help="""Will reformat the gene name field if it is split
                                  by "|" characters. The gene is in the
                                  specified 0-based column""",
                          default=None)
    opt_parser.add_option("--fpkms",
                          dest="fpkms",
                          type="string",
                          help="""Optional: Output file listing all fpkms from the input
                                  table.""",
                          default=None)
    opt_parser.add_option("--min_fpkm",
                          dest="min_fpkm",
                          type="float",
                          help="""Minimum FPKM value for a given gene. If the
                                  gene is expressed below this level in all
                                  samples, the gene is filtered from the
                                  table. DEF=%.2f""" % DEF_MIN_FPKM,
                          default=DEF_MIN_FPKM)
    opt_parser.add_option("--min_fold_fpkm",
                          dest="min_fold_fpkm",
                          type="float",
                          help="""If table will be used for differential
                                  expression analysis, then the minimum fold
                                  difference between the minimum fpkm and
                                  maximum fpkm can be used to remove genes.
                                  DEF=%.2f""" % DEF_MIN_FOLD_FPKM,
                          default=DEF_MIN_FOLD_FPKM)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--in_table")
    opt_parser.check_required("--out_table")

    in_table = open(options.in_table)
    out_table = open(options.out_table, "w")

    gene_col = options.gene_col

    fpkm_out = None
    if options.fpkms:
        fpkm_out = open(options.fpkms, "w")

    min_fpkm = options.min_fpkm
    min_fold_fpkm = options.min_fold_fpkm

    start_idx = options.start_idx

    log2_transform = options.log2_transform

    for line in in_table:
        if line.startswith("#"):
            out_table.write(line)
            continue
   
        line = formatLine(line)

        fpkms = map(float, line.split("\t")[start_idx:])

        if not check_min_fpkm(fpkms, min_fpkm):
            continue 

        fpkm_list = transform(fpkms, log2_transform)

        if fpkm_out:
            for fpkm in fpkm_list:
                fpkm_out.write(repr(fpkm) + "\n")

        if not check_min_fold_fpkm(fpkm_list, min_fold_fpkm):
            continue

        outlineList = line.split("\t")[:start_idx]
        for fpkm in fpkm_list:
            outlineList.append("%.6f" % fpkm)
  
        if not gene_col is None:
            geneName = outlineList[gene_col]
            outlineList[gene_col] = geneName.split("|")[0]

        # If all filters passed, then print line
        out_table.write("\t".join(outlineList) + "\n")

    in_table.close()
    out_table.close()

    if fpkm_out:
        fpkm_out.close()
        
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def check_min_fold_fpkm(fpkms, min_fold_fpkm):
    min_fpkm = INFINITY
    max_fpkm = -INFINITY

    for fpkm in fpkms:
        if fpkm < min_fpkm:
            min_fpkm = fpkm
        if fpkm > max_fpkm:
            max_fpkm = fpkm

    if min_fpkm == 0.0:
        return True

    fc = abs(max_fpkm/min_fpkm)

    if fc >= min_fold_fpkm:
        return True

    return False
def check_min_fpkm(fpkms, min_fpkm):
    for fpkm in fpkms:
        if fpkm >= min_fpkm:
            return True
    return False        

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def transform(fpkm_list, log2transform):
    
    
    if log2transform:
        transformed_fpkms = []
        for fpkm in fpkm_list:
            if fpkm == 0.0:
                transformed_fpkms.append(-10)
            else:
                transformed_fpkms.append(math.log(fpkm,2))
    else:
        transformed_fpkms = fpkm_list

    return transformed_fpkms
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
