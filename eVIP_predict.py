#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
#!/usr/bin/python
# mutation_impact
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
import random
import rpy2.robjects as robjects
import cmap.io.gct as gct
import cmap.io.plategrp as grp
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv

#############
# CONSTANTS #
#############
NUM_ITERATIONS = 1000
NUM_REPS = 3
LOG10_ZERO = 35.0
WT_IDX = 5
DEF_IE_COL = "x_ie_a549"
DEF_ALLELE_COL = "x_mutation_status"
DEF_IE_FILTER = 0

#predict
CONN_THRESH = 0.05
MUT_WT_THRESH = 0.05
DISTING_THRESH = 0.05
DIFF_WT_MT_RANK = 0


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
    opt_parser.add_option("--sig_info",
                          dest="sig_info",
                          type="string",
                          help="""sig info file with gene information and distil
                                  information""",
                          default=None)
    opt_parser.add_option("--gctx",
                          dest="gctx",
                          type="string",
                          help="""GCTX file of pairwise similarity. Designed to
                                  use connectivity score, but any measurement of
                                  similarity can be used (e.g., pearson
                                  correlation)""",
                          default=None)
    opt_parser.add_option("--allele_col",
                          dest="allele_col",
                          type="string",
                          help="""Column name that indicates the allele names.
                                  DEF=%s""" % DEF_ALLELE_COL,
                          default=DEF_ALLELE_COL)
    opt_parser.add_option("-o",
                          dest="output_table",
                          type="string",
                          help="""Prefix of output files of mutation impact data.
                                  Includes figures of p-value distribution.""",
                          default=None)
    opt_parser.add_option("-r",
                          dest="reference_test_file",
                          type="string",
                          help="""File explicitly indicating which comparisons
                                  to do. Assumes the file has a header and it is
                                  ignored. The first column is the reference
                                  allele and second column is test allele. If
                                  this file is not given, then the reference
                                  alleles are assumed to be WT and inferred from
                                  the allele names.""",
                          default=None)
    opt_parser.add_option("-c",
                          dest="controls_file",
                          type="string",
                          help=""".grp file containing allele names of control
                                  perturbations. If this file is given, a null
                                  will be calculated from these""",
                          default=None)
    opt_parser.add_option("-i",
                          dest="num_iterations",
                          type="int",
                          help="Number of iterations to run. DEF=%d" % NUM_ITERATIONS,
                          default=NUM_ITERATIONS)
    # opt_parser.add_option("--rep_null",
    #                       dest="rep_null_input",
    #                       type="string",
    #                       help="""Optional file containing rep null values from a previous
    #                               run. Should end in _rep_null.txt""",
    #                      default=None)
    opt_parser.add_option("--conn_null",
                          dest="conn_null_input",
                          type="string",
                          help="""Optional file containing connectvity null values from a previous
                                  run. Should end in _conn_null.txt""",
                          default=None)
    opt_parser.add_option("--ie_col",
                          dest="ie_col",
                          type="string",
                          help="""Name of the column with infection efficiency
                                  information. DEF=%s""" % DEF_IE_COL,
                          default=DEF_IE_COL)
    opt_parser.add_option("--ie_filter",
                          dest="ie_filter",
                          type="float",
                          help="""Threshold for infection efficiency. Any wildtype
                                  or mutant alleles having an ie below this
                                  threshold, will be removed""",
                          default=DEF_IE_FILTER)
    opt_parser.add_option("--num_reps",
                          dest="num_reps",
                          type="int",
                          help="""Number of replicates expected for each allele.
                                  DEF=%d""" % NUM_REPS,
                          default=NUM_REPS)
    opt_parser.add_option("--cell_id",
                          dest="cell_id",
                          type="string",
                          help="""Optional: Will only look at signatures from this cell
                                  line. Helps to filter sig_info file.""",
                          default=None)
    opt_parser.add_option("--plate_id",
                          dest="plate_id",
                          type="string",
                          help="""Optional: Will only look at signatures from
                                  this plate.""",
                          default=None)

    #predict
    opt_parser.add_option("--conn_thresh",
                          dest="conn_thresh",
                          type="float",
                          help="P-value threshould for connectivity vs null.",
                          default=CONN_THRESH)
    opt_parser.add_option("--mut_wt_rep_thresh",
                          dest="mut_wt_thresh",
                          type="float",
                          help="""P-value threshould for comparison of WT and mut
                                  robustness""",
                          default=MUT_WT_THRESH)
    opt_parser.add_option("--mut_wt_rep_rank_diff",
                          dest="mut_wt_rep_diff",
                          type="float",
                          help="""The minimum difference in median rankpoint
                                  between WT and mut to consider a difference.
                                  DEF=%d""" % DIFF_WT_MT_RANK,
                          default=DIFF_WT_MT_RANK)
    opt_parser.add_option("--disting_thresh",
                          dest="disting_thresh",
                          type="float",
                          help="""P-value threshould that tests if mut and wt reps
                                  are indistinguishable from each other""",
                          default=DISTING_THRESH)
    opt_parser.add_option("--use_c_pval",
                          dest="use_c_pval",
                          action="store_true",
                          help="Will use corrected p-value instead of raw p-val",
                          default=False)
    # opt_parser.add_option("--conn_null_med",
    #                   dest="conn_null_med",
    #                   type="float",
    #                   help="Median of null connectivity distribution.",
    #                   default=CONN_NULL_MED)

    (options, args) = opt_parser.parse_args()

    #validate the command line arguments
    opt_parser.check_required("--sig_info")
    opt_parser.check_required("--gctx")
    opt_parser.check_required("-c")
    opt_parser.check_required("-o")
    #validate from predict
    opt_parser.check_required("--mut_wt_rep_thresh")
    opt_parser.check_required("--disting_thresh")

    sig_info_file = open(options.sig_info)
    output = open(options.output_table, "w")


    # Output distribution files
    controls = grp.read_grp(options.controls_file)

    reference_test_filename = options.reference_test_file
    ref2test_allele = None
    if reference_test_filename:
        ref2test_allele = parseRefTestFile(reference_test_filename)

    allele_col = options.allele_col

    this_gctx = gct.GCT(options.gctx)
    this_gctx.read()

    num_iterations = options.num_iterations
    num_reps = options.num_reps

    ie_col = options.ie_col
    ie_filter = options.ie_filter

    #rep_null_input = options.rep_null_input
    conn_null_input = options.conn_null_input

    #options from predict
    c_thresh = options.conn_thresh
    mut_wt_thresh = options.mut_wt_thresh
    disting_thresh = options.disting_thresh
    mut_wt_rep_diff = options.mut_wt_rep_diff
    use_c_pval = options.use_c_pval
    #conn_null_med = options.conn_null_med


    if conn_null_input:
        conn_nulls_from_input_str = grp.read_grp(conn_null_input)
        conn_nulls_from_input = map(float, conn_nulls_from_input_str)

    (allele2distil_id,
     allele2WT,
     allele2gene,
     allele2cell_id,
     WT_alleles) = parse_sig_info(sig_info_file,
                                  ref2test_allele,
                                  allele_col,
                                  ie_col, ie_filter,
                                  options.cell_id,
                                  options.plate_id)


    clean_controls = []
    for this_control in controls:
        if this_control in allele2distil_id:
            clean_controls.append(this_control)

    if not conn_null_input:
        replicate_null_dist, connectivity_null_dist = getNullDist(this_gctx,
                                                              allele2distil_id,
                                                              clean_controls,
                                                              num_iterations,
                                                              num_reps)

    if conn_null_input:
        connectivity_null_dist = conn_nulls_from_input


    if not conn_null_input:
        conn_null_dist_out = open(options.output_table+ "_conn_null.txt", "w")
        for x in connectivity_null_dist:
            conn_null_dist_out.write("%f\n" % x)
        conn_null_dist_out.close()


    WT_dict = buildWT_dict(this_gctx, allele2distil_id, WT_alleles, num_reps)

    column_headers = ["gene",
                      "mut",
                      "mut_rep",
                      "wt_rep",
                      "mut_wt_connectivity",
                      "wt",
                      "cell_line",
                      "mut_wt_rep_pval",
                      "mut_wt_conn_null_pval",
                      "wt_mut_rep_vs_wt_mut_conn_pval",
                      "mut_wt_rep_c_pval",
                      "mut_wt_conn_null_c_pval",
                      "wt_mut_rep_vs_wt_mut_conn_c_pval",
                      "prediction"]

    file_writer = csv.DictWriter(output, delimiter="\t",fieldnames=column_headers)
    file_writer.writeheader()

    mut_wt_rep_pvals = []
    mut_wt_conn_pvals = []
    mut_wt_rep_vs_wt_mut_conn_pvals = []
    outlines = []
    predictions =[]

    # Build comparison
    for allele in allele2WT:

        # Don't calculate for the WT allele
        if allele == allele2WT[allele]:
            continue

        mut_rankpt, mut_rankpt_dist = getSelfConnectivity(this_gctx,
                                                          allele2distil_id[allele],
                                                          num_reps)

        mut_wt_conn_rankpt, mut_wt_conn_dist = getConnectivity(this_gctx,
                                                               allele2distil_id[allele],
                                                               allele2distil_id[allele2WT[allele]],
                                                               num_reps)

        conn_pval = getPairwiseComparisons(mut_wt_conn_dist,
                                           connectivity_null_dist)
        mut_wt_conn_pvals.append(conn_pval)

        mut_wt_rep_pval = getPairwiseComparisons(mut_rankpt_dist,
                                                 WT_dict[allele2WT[allele]]["wt_rep_dist"])
        mut_wt_rep_pvals.append(mut_wt_rep_pval)

        wt_mut_rep_vs_wt_mut_conn_pval = getKruskal(WT_dict[allele2WT[allele]]["wt_rep_dist"],
                                                    mut_rankpt_dist,
                                                    mut_wt_conn_dist)
        mut_wt_rep_vs_wt_mut_conn_pvals.append(wt_mut_rep_vs_wt_mut_conn_pval)

        # Calculate corrected pvalues
        mut_wt_rep_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_rep_pvals), "BH")
        mut_wt_conn_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_conn_pvals), "BH")
        mut_wt_rep_vs_wt_mut_conn_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_rep_vs_wt_mut_conn_pvals), "BH")

        out_elems = [allele2gene[allele],
                     allele,
                     "%f" % mut_rankpt,
                     "%f" % WT_dict[allele2WT[allele]]["wt_rep"],
                     "%f" % mut_wt_conn_rankpt,
                     allele2WT[allele],
                     allele2cell_id[allele],
                     "%f" % mut_wt_rep_pval,
                     "%f" % conn_pval,
                     "%f" % wt_mut_rep_vs_wt_mut_conn_pval]


        outline = "\t".join(out_elems)
        outlines.append(outline)

        if use_c_pval:
            prediction = get_prediction_6(float(WT_dict[allele2WT[allele]]["wt_rep"]),
                                          float(mut_rankpt),
                                          float(mut_wt_rep_c_pvals[i]),
                                          float(mut_wt_conn_rankpt),
                                          float(mut_wt_conn_c_pvals[i]),
                                          float(mut_wt_rep_vs_wt_mut_conn_c_pvals[i]),
                                          mut_wt_thresh,
                                          mut_wt_rep_diff,
                                          c_thresh,
                                          disting_thresh)
            predictions.append(prediction)
        else:
            prediction = get_prediction_6(float(WT_dict[allele2WT[allele]]["wt_rep"]),
                                          float(mut_rankpt),
                                          float(mut_wt_rep_pval),
                                          float(mut_wt_conn_rankpt),
                                          float(conn_pval),
                                          float(wt_mut_rep_vs_wt_mut_conn_pval),
                                          mut_wt_thresh,
                                          mut_wt_rep_diff,
                                          c_thresh,
                                          disting_thresh)

            predictions.append(prediction)


    # Write to file
    #how many outlines?
    num_lines = len(outlines)
    count = 0
    #for each outline
    for i in range(num_lines):
        this_outline = outlines[i]
        # Getting wt c_pval
        this_outlist = outlines[i].split("\t")
        this_wt = this_outlist[WT_IDX]
        #wt_idx = wt_ordered.index(this_wt)

        this_outline += "\t%f\t" % mut_wt_rep_c_pvals[i]
        this_outline += "%f\t" % mut_wt_conn_c_pvals[i]
        this_outline += "%f\t" % mut_wt_rep_vs_wt_mut_conn_c_pvals[i]
        this_outline += "%s\t" % predictions[i]
        this_outline += "\n"

        output.write(this_outline)

def eVIP_run_main(sig_info=None, o=None, c=None, r=None, gctx=None, conn_thresh=None, conn_null=None,
                  allele_col=None,
                  ie_col=None, ie_filter=None, cell_id=None, plate_id=None,
                  i=None, num_reps=None, mut_wt_rep_thresh=None,
                  disting_thresh=None, mut_wt_rep_rank_diff=None, use_c_pval=None):

    if sig_info == None:
        raise Exception("Missing sig_info file in function call")
    if o == None:
        raise Exception("Missing input file in function call")
    if c == None:
        raise Exception("Missing input file in function call")
    if r == None:
        raise Exception("Missing input file in function call")
    if gctx == None:
        raise Exception("Missing input file in function call")

    # setting defaults (double check if numbers are right)
    ie_filter = float(ie_filter) if ie_filter != None else float(0.0)
    conn_thresh = float(conn_thresh) if conn_thresh != None else float(0.05)
    allele_col = str(allele_col) if allele_col != None else str(x_mutation_status)
    ie_col = str(ie_col) if ie_col != None else str(x_ie_a549)
    # do i need the "" here?
    i = int(i) if i != None else int(1000)
    num_reps = int(num_reps) if num_reps != None else int(3)
    mut_wt_rep_thresh = float(mut_wt_rep_thresh) if mut_wt_rep_thresh != None else float(0.05)
    disting_thresh = float(disting_thresh) if disting_thresh != None else float(0.05)
    mut_wt_rep_rank_diff = float(mut_wt_rep_rank_diff) if mut_wt_rep_rank_diff != None else float(0)
    #conn_null_med = float(conn_null_med) if conn_null_med != None else float(0.05)

    sig_info_file = open(sig_info)
    output = open(o+".txt", "w")

    # Output distribution files
    controls = grp.read_grp(c)

    reference_test_filename = r
    ref2test_allele = None
    if reference_test_filename:
        ref2test_allele = parseRefTestFile(reference_test_filename)

    this_gctx = gct.GCT(gctx)
    this_gctx.read()

    num_iterations = i
    c_thresh = conn_thresh
    num_reps = int(num_reps)
    mut_wt_thresh = mut_wt_rep_thresh
    mut_wt_rep_diff = mut_wt_rep_rank_diff


    if conn_null:
        conn_nulls_from_input_str = grp.read_grp(conn_null)
        conn_nulls_from_input = map(float, conn_nulls_from_input_str)


    (allele2distil_id,
     allele2WT,
     allele2gene,
     allele2cell_id,
     WT_alleles) = parse_sig_info(sig_info_file,
                                  ref2test_allele,
                                  allele_col,
                                  ie_col,
                                  ie_filter,
                                  cell_id,
                                  plate_id)

    clean_controls = []

    for this_control in controls:
        if this_control in allele2distil_id:
            clean_controls.append(this_control)

    if not conn_null:
        replicate_null_dist, connectivity_null_dist = getNullDist(this_gctx,
                                                                  allele2distil_id,
                                                                  clean_controls,
                                                                  num_iterations,
                                                                  num_reps)

    if conn_null:
        connectivity_null_dist = conn_nulls_from_input

    if not conn_null:
        conn_null_dist_out = open(o + "_conn_null.txt", "w")
        for x in connectivity_null_dist:
            conn_null_dist_out.write("%f\n" % x)
        conn_null_dist_out.close()

    WT_dict = buildWT_dict(this_gctx, allele2distil_id, WT_alleles, num_reps)

    column_headers = ["gene",
                      "mut",
                      "mut_rep",
                      "wt_rep",
                      "mut_wt_connectivity",
                      "wt",
                      "cell_line",
                      "mut_wt_rep_pval",
                      "mut_wt_conn_null_pval",
                      "wt_mut_rep_vs_wt_mut_conn_pval",
                      "mut_wt_rep_c_pval",
                      "mut_wt_conn_null_c_pval",
                      "wt_mut_rep_vs_wt_mut_conn_c_pval",
                      "prediction"]

    file_writer = csv.DictWriter(output, delimiter="\t", fieldnames=column_headers)
    file_writer.writeheader()

    mut_wt_rep_pvals = []
    mut_wt_conn_pvals = []
    mut_wt_rep_vs_wt_mut_conn_pvals = []
    outlines = []
    predictions = []

    # Build comparison
    for allele in allele2WT:

        # Don't calculate for the WT allele
        if allele == allele2WT[allele]:
            continue

        mut_rankpt, mut_rankpt_dist = getSelfConnectivity(this_gctx,
                                                          allele2distil_id[allele],
                                                          num_reps)

        mut_wt_conn_rankpt, mut_wt_conn_dist = getConnectivity(this_gctx,
                                                               allele2distil_id[allele],
                                                               allele2distil_id[allele2WT[allele]],
                                                               num_reps)

        conn_pval = getPairwiseComparisons(mut_wt_conn_dist,
                                           connectivity_null_dist)

        mut_wt_conn_pvals.append(conn_pval)

        mut_wt_rep_pval = getPairwiseComparisons(mut_rankpt_dist,
                                                 WT_dict[allele2WT[allele]]["wt_rep_dist"])
        mut_wt_rep_pvals.append(mut_wt_rep_pval)

        wt_mut_rep_vs_wt_mut_conn_pval = getKruskal(WT_dict[allele2WT[allele]]["wt_rep_dist"],
                                                    mut_rankpt_dist,
                                                    mut_wt_conn_dist)
        mut_wt_rep_vs_wt_mut_conn_pvals.append(wt_mut_rep_vs_wt_mut_conn_pval)

        # Calculate corrected pvalues
        mut_wt_rep_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_rep_pvals), "BH")
        mut_wt_conn_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_conn_pvals), "BH")
        mut_wt_rep_vs_wt_mut_conn_c_pvals = robjects.r['p.adjust'](
            robjects.FloatVector(mut_wt_rep_vs_wt_mut_conn_pvals), "BH")

        out_elems = [allele2gene[allele],
                     allele,
                     "%f" % mut_rankpt,
                     "%f" % WT_dict[allele2WT[allele]]["wt_rep"],
                     "%f" % mut_wt_conn_rankpt,
                     allele2WT[allele],
                     allele2cell_id[allele],
                     "%f" % mut_wt_rep_pval,
                     "%f" % conn_pval,
                     "%f" % wt_mut_rep_vs_wt_mut_conn_pval]

        outline = "\t".join(out_elems)
        outlines.append(outline)

        if use_c_pval:
            prediction = get_prediction_6(float(WT_dict[allele2WT[allele]]["wt_rep"]),
                                          float(mut_rankpt),
                                          float(mut_wt_rep_c_pvals[i]),
                                          float(mut_wt_conn_rankpt),
                                          float(mut_wt_conn_c_pvals[i]),
                                          float(mut_wt_rep_vs_wt_mut_conn_c_pvals[i]),
                                          mut_wt_thresh,
                                          mut_wt_rep_diff,
                                          c_thresh,
                                          disting_thresh)
            predictions.append(prediction)
        else:
            prediction = get_prediction_6(float(WT_dict[allele2WT[allele]]["wt_rep"]),
                                          float(mut_rankpt),
                                          float(mut_wt_rep_pval),
                                          float(mut_wt_conn_rankpt),
                                          float(conn_pval),
                                          float(wt_mut_rep_vs_wt_mut_conn_pval),
                                          mut_wt_thresh,
                                          mut_wt_rep_diff,
                                          c_thresh,
                                          disting_thresh)

            predictions.append(prediction)

    # Write to file
    # how many outlines?
    num_lines = len(outlines)
    count = 0
    # for each outline
    for i in range(num_lines):
        this_outline = outlines[i]
        # Getting wt c_pval
        this_outlist = outlines[i].split("\t")
        this_wt = this_outlist[WT_IDX]
        # wt_idx = wt_ordered.index(this_wt)

        this_outline += "\t%f\t" % mut_wt_rep_c_pvals[i]
        this_outline += "%f\t" % mut_wt_conn_c_pvals[i]
        this_outline += "%f\t" % mut_wt_rep_vs_wt_mut_conn_c_pvals[i]
        this_outline += "%s\t" % predictions[i]
        this_outline += "\n"

        output.write(this_outline)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def buildWT_dict(this_gctx, allele2distil_id, WT_alleles, num_reps):
    """
    {WT_allele:{"wt_rep": med_wt_rep,
                "wt_rep_dist":[]
                "wt_rep_pval": p_val vs null}
    """
    WT_dict = {}

    for allele in WT_alleles:
        WT_dict[allele] = {}
        wt_rep_rankpt, rep_rankpts = getSelfConnectivity(this_gctx,
                                                        allele2distil_id[allele],
                                                        num_reps)
        WT_dict[allele]["wt_rep"] = wt_rep_rankpt
        WT_dict[allele]["wt_rep_dist"] = rep_rankpts

    return WT_dict

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getKruskal(wt_rankpt_dist, mut_rankpt_dist, mut_wt_conn_dist):
    return robjects.r["kruskal.test"](robjects.ListVector({'a':robjects.FloatVector(wt_rankpt_dist),
                                                           'b':robjects.FloatVector(mut_rankpt_dist),
                                                           'c':robjects.FloatVector(mut_wt_conn_dist)}))[2][0]


def getSelfConnectivity(this_gctx, distil_ids, num_reps):
    """
    returns
    median rankpoint
    Distribution of values
    """
    row_medians = []
    for i in range(num_reps):
        rank_pts = []
        for j in range(num_reps):
            if i == j:
                continue
            rank_pts.append(float(this_gctx.frame[distil_ids[i]]
                                                 [distil_ids[j]]))

        row_medians.append(numpy.percentile(rank_pts, 50))

    return numpy.percentile(row_medians, 50), row_medians


#   rep_rankpts = []
#   for i in range(num_reps):
#       for j in range(i+1, num_reps):
#           rep_rankpts.append(float(this_gctx.frame[distil_ids[i]]
#                                                   [distil_ids[j]]))
#
#   return numpy.percentile(rep_rankpts, 50), rep_rankpts

def getConnectivity(this_gctx, distil_ids1, distil_ids2, num_reps):
    """
    Returns
    median connectivity
    Distribution of values
    """
    row_column_medians = []

    col_rankpts = [[] for n in range(num_reps)]
    for i in range(num_reps):
        row_rankpts = []
        for j in range(num_reps):
            val = float(this_gctx.frame[distil_ids1[i]][distil_ids2[j]])
            row_rankpts.append(val)
            col_rankpts[j].append(val)

        row_column_medians.append(numpy.percentile(row_rankpts, 50))

    for n in range(num_reps):
        row_column_medians.append(numpy.percentile(col_rankpts[n], 50))

    return numpy.percentile(row_column_medians, 50), row_column_medians


def getLog(vals):
    out_vals = []
    for val in vals:
        if val == 0:
            out_vals.append(LOG10_ZERO)
            continue
        out_vals.append(-math.log10(val))

    return out_vals

def getPairwiseComparisons(dist1, dist2):

    return robjects.r["wilcox.test"](robjects.FloatVector(dist1),
                                     robjects.FloatVector(dist2))[2][0]

def getNullDist(this_gctx, allele2distil_id, controls,
                num_iterations, num_reps):
    """
    Returns
    replicate_null_dist,
    connectivity_null_dist
    """
    rep_null_dist = []
    connectivity_null_dist = []

    allele_list = allele2distil_id.keys()



    c = 0
    while c < num_iterations:
        random_control = random.choice(controls)
        random_allele = random.choice(allele_list)

        if random_control == random_allele:
            continue

        control_distil_ids_set = set([random.choice(allele2distil_id[random_control])])
        while len(control_distil_ids_set) < num_reps:
            this_control = random.choice(controls)
            this_distil = random.choice(allele2distil_id[this_control])
            if this_distil in control_distil_ids_set:
                continue
            control_distil_ids_set.add(this_distil)

        control_distil_ids = list(control_distil_ids_set)


        # Introspect similarity is median of row similarities
        rank_pts = []
        for i in range(1, num_reps):
            rank_pts.append(float(this_gctx.frame[control_distil_ids[0]][control_distil_ids[i]]))


        rep_null_dist.append(numpy.percentile(rank_pts, 50))

        # Get connectivity null
        rank_pts = []
        for i in range(num_reps):
            rank_pts.append(float(this_gctx.frame[allele2distil_id[random_control][0]]
                                                 [allele2distil_id[random_allele][i]]))


        connectivity_null_dist.append(numpy.percentile(rank_pts, 50))

        c += 1

    return rep_null_dist, connectivity_null_dist


def hasLowIE(ie_string, ie_thresh):
    ie_elems = map(float, ie_string.split("|"))
    for ie in ie_elems:
        if ie < ie_thresh:
            return True
    return False


def parseRefTestFile(reference_test_filename):

    ref2test = {}

    ref_test_file = open(reference_test_filename)
    ctr = 1
    for line in ref_test_file:
        if ctr == 1:
            ctr += 1
            continue

        line = formatLine(line)
        lineList = line.split("\t")

        updateDictOfLists(ref2test, lineList[0], lineList[1])

    return ref2test

def parse_sig_info(sig_info_file, ref2test_allele, allele_col, ie_col, ie_filter = None, cell_id = None, plate_id = None):
    """
    Returns:
    allele2distil_id = {}
    allele2WT = {}
    allele2gene,
    allele2cell_id
    WT_alleles = []
    """
    allele2distil_id = {}
    gene2alleles = {}
    allele2WT = {}

    allele2gene = {}
    allele2cell_id = {}

    ie_col_idxs = []

    passed_alleles = set([])

    for line in sig_info_file:
        line = formatLine(line)
        lineList = line.split("\t")

        if "distil_id" in line:
            sig_id_idx = lineList.index("sig_id")
            distil_idx = lineList.index("distil_id")
            gene_idx = lineList.index("pert_mfc_desc")
            allele_idx = lineList.index(allele_col)
            cell_idx = lineList.index("cell_id")

            ie_cols = ie_col.split(",")
            for ie_col in ie_cols:
                ie_col_idxs.append(lineList.index(ie_col))
            continue

        if cell_id:
            if cell_id != lineList[cell_idx]:
                continue

        if plate_id:
            if plate_id not in lineList[sig_id_idx]:
                continue

        if ie_filter:
            for ie_idx in ie_col_idxs:
                low_ie_flag = hasLowIE(lineList[ie_idx], ie_filter)
                if low_ie_flag:
                    break
            if low_ie_flag:
                continue

        x_mutation_status = lineList[allele_idx]
        gene = lineList[gene_idx]

        if x_mutation_status in allele2distil_id:
            print "Warning: Duplicate allele present %s" % x_mutation_status
            continue

        distil_ids = lineList[distil_idx].split("|")

        for distil_id in distil_ids:
            updateDictOfLists(allele2distil_id, x_mutation_status, distil_id)

        updateDictOfSets(gene2alleles, gene, x_mutation_status)

        allele2gene[x_mutation_status] = gene
        allele2cell_id[x_mutation_status] = lineList[cell_idx]
        passed_alleles.add(x_mutation_status)


    if ref2test_allele:
        WT_alleles = set(ref2test_allele.keys())
        WT_alleles = WT_alleles & passed_alleles

        for ref in ref2test_allele:
            if ref not in passed_alleles:
                continue
            for test in ref2test_allele[ref]:
                if test not in passed_alleles:
                    continue
                allele2WT[test] = ref
    else: # Need to infer WT reference
        # Now find WT orf
        WT_alleles = set([])

        for gene in gene2alleles:
            # Find WT allele
            this_WT = ""
            if gene + "_WT.c" in gene2alleles[gene]:
                this_WT = gene + "_WT.c"
            elif gene + "_WT.c.2" in gene2alleles[gene]:
                this_WT = gene + "_WT.c.2"
            elif gene + "_WT.o" in gene2alleles[gene]:
                this_WT = gene + "_WT.o"
            elif gene + "_WT" in gene2alleles[gene]:
                this_WT = gene + "_WT"
            elif gene + "_WT.1" in gene2alleles[gene]:
                this_WT = gene + "_WT.1"
            elif gene + "_WT.2" in gene2alleles[gene]:
                this_WT = gene + "_WT.2"
            else:
                print "Warning: Can't find WT orf for %s" % gene
                continue

            WT_alleles.add(this_WT)

            # Match up alleles to WT
            for allele in gene2alleles[gene]:
                allele2WT[allele] = this_WT

    return allele2distil_id, allele2WT, allele2gene, allele2cell_id, list(WT_alleles)

def updateDictOfLists(d, key, item):
    """
    """
    try:
        d[key].append(item)
    except KeyError:
        d[key] = [item]

def updateDictOfSets(d, key, item):
    """
    Similar functionality as updateDictOfLists
    """
    try:
        d[key].add(item)
    except KeyError:
        d[key] = set([item])
def get_prediction_6(wt_rep, mut_rep, mut_wt_rep_pval,
                     mut_wt_conn, mut_wt_conn_pval, disting_pval,
                     mut_wt_thresh, mut_wt_rep_diff, c_thresh, disting_thresh):

    if disting_pval < disting_thresh:
        if max_diff(wt_rep, mut_rep, mut_wt_conn) < mut_wt_rep_diff:
            if mut_wt_conn_pval < c_thresh:
                return "Neutral"
            else:
                return "Error"

#       if mut_wt_conn_pval < c_thresh:
#           if mut_wt_conn < conn_null_med:
#               return "DOM-NEG"

        if mut_wt_rep_pval < mut_wt_thresh:
            if wt_rep < mut_rep:
                if mut_rep - wt_rep >= mut_wt_rep_diff:
                    return "GOF"
                else:
                    return "COF"
            elif wt_rep > mut_rep:
                if wt_rep - mut_rep >= mut_wt_rep_diff:
                    return "LOF"
                else:
                    return "COF"
            else:
                return "COF"
        else:
            return "COF"

    if mut_wt_conn_pval < c_thresh:
        return "Neutral"

    return "NI"

def max_diff(wt_rep, mut_rep, mut_wt_conn):
    max_diff = abs(wt_rep - mut_rep)

    wt_conn_diff = abs(wt_rep - mut_wt_conn)
    if wt_conn_diff > max_diff:
        max_diff = wt_conn_diff

    mut_conn_diff = abs(mut_rep - mut_wt_conn)
    if mut_conn_diff > max_diff:
        max_diff = mut_conn_diff

    return max_diff


#################
# END FUNCTIONS #
#################
if __name__ == "__main__": main()
