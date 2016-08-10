#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
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

#############
# CONSTANTS #
#############
NUM_ITERATIONS = 1000
NUM_REPS = 3
IE_FILTER = 0

#LOG10_ZERO = 10.0
LOG10_ZERO = 35.0

WT_IDX = 5

DEF_IE_COL = "x_ie_a549"
DEF_ALLELE_COL = "x_mutation_status"

CONN_THRESH = 0.05
MUT_WT_REP_THRESH = 0.05
DISTING_THRESH = 0.05
DIFF_WT_MT_RANK = 5.0

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
                          dest="output_file_prefix",
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
#   opt_parser.add_option("--null",
#                         dest="out_null_file",
#                         type="string",
#                         help="""File of pairwise comparisons between alleles
#                                 from different genes""",
#                         default=None)
    opt_parser.add_option("-i",
                          dest="num_iterations",
                          type="int",
                          help="Number of iterations to run. DEF=%d" % NUM_ITERATIONS,
                          default=NUM_ITERATIONS)
#    opt_parser.add_option("--rep_null",
#                          dest="rep_null_input",
#                          type="string",
#                          help="""Optional file containing rep null values from a previous
#                                  run. Should end in _rep_null.txt""",
#                          default=None)
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
                                  threshold, will be removed. DEF=%d""" % IE_FILTER,
                          default=IE_FILTER)
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
                                  line. Helps to filter sig_info file. DEF=%s""",
                          default=None)
    opt_parser.add_option("--plate_id",
                          dest="plate_id",
                          type="string",
                          help="""Optional: Will only look at signatures from
                                  this plate.""",
                          default=None)

    #prediction options

    opt_parser.add_option("--use_c_pval",
                          dest="use_c_pval",
                          action="store_true",
                          help="Will use corrected p-value instead of raw p-val",
                          default=False)
    opt_parser.add_option("--conn_thresh",
                          dest="conn_thresh",
                          type="float",
                          help="P-value threshould for connectivity vs null. DEF=%d" % CONN_THRESH,
                          default=CONN_THRESH)
    opt_parser.add_option("--mut_wt_rep_thresh",
                          dest="mut_wt_thresh",
                          type="float",
                          help="""P-value threshould for comparison of WT and mut
                                  robustness DEF=%d""" % MUT_WT_REP_THRESH,
                          default=MUT_WT_REP_THRESH)
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
                                  are indistinguishable from each other. DEF=%d""" % DISTING_THRESH,
                          default=DISTING_THRESH)

    (options, args) = opt_parser.parse_args()

    # validate the command line arguments
    opt_parser.check_required("--sig_info")
    opt_parser.check_required("--gctx")
#    opt_parser.check_required("--null")
    opt_parser.check_required("-c")
    opt_parser.check_required("-o")

    sig_info_file = open(options.sig_info)
#    out_null_file = open(options.out_null_file, "w")
    output_file_prefix = open(options.output_file_prefix + ".txt", "w")

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

#    rep_null_input = options.rep_null_input
    conn_null_input = options.conn_null_input


#predictions
    use_c_pval = options.use_c_pval
    c_thresh = options.conn_thresh
    mut_wt_thresh = options.mut_wt_thresh
    mut_wt_rep_diff = options.mut_wt_rep_diff
    disting_thresh = options.disting_thresh



#    if rep_null_input:
#        rep_nulls_from_input_str = grp.read_grp(rep_null_input)
#        rep_nulls_from_input = map(float, rep_nulls_from_input_str)
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

#    if (not rep_null_input) or (not conn_null_input):
    replicate_null_dist, connectivity_null_dist = getNullDist(this_gctx,
                                                              allele2distil_id,
                                                              clean_controls,
                                                              num_iterations,
                                                              num_reps)

#    if rep_null_input:
#        replicate_null_dist = rep_nulls_from_input
    if conn_null_input:
        connectivity_null_dist = conn_nulls_from_input

    # Print out percentiles of each null distribution
    # print "Replicate null percentiles"
    #print "2.5,5,10,50,90,95,97.5"
    rep_precentiles =  numpy.percentile(replicate_null_dist, [2.5,5,10,50,90,95,97.5])
    #if not rep_null_input:
    rep_null_distribution_out = open(options.output_file_prefix + "_rep_null.txt", "w")
    for x in replicate_null_dist:
        rep_null_distribution_out.write("%f\n" % x)
    rep_null_distribution_out.close()
    #print rep_precentiles

    #print "Connectivity null percentiles"
    #print "2.5,5,10,50,90,95,97.5"
    conn_percentiles = numpy.percentile(connectivity_null_dist, [2.5,5,10,50,90,95,97.5])
    if not conn_null_input:
        conn_null_dist_out = open(options.output_file_prefix + "_conn_null.txt", "w")
        for x in connectivity_null_dist:
            conn_null_dist_out.write("%f\n" % x)
        conn_null_dist_out.close()
    #print conn_percentiles

    # WT null data
    # {WT_allele:{"wt_rep": med_wt_rep,
    #             "wt_rep_dist":[]
    #             "wt_rep_pval": p_val vs null}
    WT_dict, wt_rep_pvals, wt_ordered = buildWT_dict(this_gctx, allele2distil_id, WT_alleles, replicate_null_dist,
                           num_reps)

    # Print header to output file
    output_file_prefix.write("gene\tmut\tmut_rep\twt_rep\tmut_wt_connectivity\t")
    output_file_prefix.write("wt\tcell_line\tmut_wt_rep_pval\tmut_wt_conn_null_pval\t")
    output_file_prefix.write("wt_mut_rep_vs_wt_mut_conn_pval\tmut_wt_rep_c_pval\t")
    output_file_prefix.write("mut_wt_conn_null_c_pval\twt_mut_rep_vs_wt_mut_conn_c_pval\t")
    output_file_prefix.write("prediction\n")


    mut_rep_pvals = []
    mut_wt_rep_pvals = []
    mut_wt_conn_pvals = []
    mut_wt_rep_vs_wt_mut_conn_pvals = []

    outlines = []

    # Build comparison
    for allele in allele2WT:

        # Don't calculate for the WT allele
        if allele == allele2WT[allele]:
            continue

        mut_rankpt, mut_rankpt_dist = getSelfConnectivity(this_gctx,
                                                          allele2distil_id[allele],
                                                          num_reps)

        self_pval = getPairwiseComparisons(mut_rankpt_dist,
                                           replicate_null_dist)
        mut_rep_pvals.append(self_pval)

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

#        wt_mut_rep_dist = WT_dict[allele2WT[allele]]["wt_rep_dist"] + mut_rankpt_dist
#        wt_mut_rep_dist = WT_dict[allele2WT[allele]]["wt_rep_dist"]
        # Through visualization, it makes more sense to do a Kruskal-Wallis test
        # to deterimine if the three categories are signficantly different.

#        wt_mut_rep_vs_wt_mut_conn_pval = getPairwiseComparisons(wt_mut_rep_dist, mut_wt_conn_dist)
        wt_mut_rep_vs_wt_mut_conn_pval = getKruskal(WT_dict[allele2WT[allele]]["wt_rep_dist"],
                                                    mut_rankpt_dist,
                                                    mut_wt_conn_dist)
        mut_wt_rep_vs_wt_mut_conn_pvals.append(wt_mut_rep_vs_wt_mut_conn_pval)

        out_elems = [allele2gene[allele],
                     allele,
                     "%f" % mut_rankpt,
                     "%f" % WT_dict[allele2WT[allele]]["wt_rep"],
                     "%f" % mut_wt_conn_rankpt,
                     allele2WT[allele],
                     allele2cell_id[allele],
    #                 "%f" % self_pval,
    #                 "%f" % WT_dict[allele2WT[allele]]["wt_rep_pval"],
                     "%f" % mut_wt_rep_pval,
                     "%f" % conn_pval,
                     "%f" % wt_mut_rep_vs_wt_mut_conn_pval]
        outline = "\t".join(out_elems)

        outlines.append(outline)

    # Calculate corrected pvalues
    mut_rep_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_rep_pvals), "BH")
    wt_rep_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(wt_rep_pvals), "BH")
    mut_wt_rep_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_rep_pvals), "BH")
    mut_wt_conn_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_conn_pvals), "BH")
    mut_wt_rep_vs_wt_mut_conn_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_rep_vs_wt_mut_conn_pvals), "BH")

    # Write to file
    num_lines = len(outlines)
    for i in range(num_lines):
        this_outline = outlines[i]

    #    this_outline += "\t%f\t" % mut_rep_c_pvals[i]

        # Getting wt c_pval
        this_outlist = outlines[i].split("\t")
        this_wt = this_outlist[WT_IDX]
        wt_idx = wt_ordered.index(this_wt)
    #    this_outline += "%f\t" % wt_rep_c_pvals[wt_idx]

        this_outline += "\t%f\t" % mut_wt_rep_c_pvals[i]
        this_outline += "%f\t" % mut_wt_conn_c_pvals[i]
        this_outline += "%f\n" % mut_wt_rep_vs_wt_mut_conn_c_pvals[i]

        output_file_prefix.write(this_outline)

    # Print out distribution files

    # self_rep
    fig1 = plt.figure()
    plt.hist(numpy.array(getLog(wt_rep_pvals)),
             bins = 50,
#             range = [0,4],
             histtype='stepfilled',
             normed=True,
             color='b',
             alpha=0.25,
             label="wt_rep")
    plt.hist(numpy.array(getLog(mut_rep_pvals)),
             bins = 50,
#             range = [0,4],
             histtype='stepfilled',
             normed=True,
             color='r',
             alpha=0.25,
             label="mut_rep")

    plt.xlabel("Replicate vs null -log10(p-val)")
    plt.ylabel("Relative frequency")
    plt.legend()

    fig1.savefig("%s_self_rep_pval_dist.png" % options.output_file_prefix)

    fig2 = plt.figure()

    plt.hist(numpy.array(getLog(mut_wt_rep_pvals)),
             bins = 50,
#             range = [0,2],
             histtype='stepfilled',
             normed=True,
             color='b',
             alpha=0.25,
             label="mut_vs_wt_rep")

    plt.xlabel("wt rep vs mut rep -log10(p-val)")
    plt.ylabel("Relative frequency")
    plt.legend()

    fig2.savefig("%s_mut_wt_rep_pval_dist.png" % options.output_file_prefix)


    fig3 = plt.figure()

    plt.hist(numpy.array(getLog(mut_wt_conn_pvals)),
             bins = 50,
#             range = [0,8],
             histtype='stepfilled',
             normed=True,
             color='b',
             alpha=0.25,
             label="mut_wt_conn")

    plt.xlabel("wt mut conn vs null -log10(p-val)")
    plt.ylabel("Relative frequency")
    plt.legend()

    fig3.savefig("%s_mut_wt_conn_pval_dist.png" % options.output_file_prefix)

    fig4 = plt.figure()

    plt.hist(numpy.array(getLog(mut_wt_rep_vs_wt_mut_conn_pvals)),
             bins = 50,
#             range = [0,4],
             histtype='stepfilled',
             normed=True,
             color='b',
             alpha=0.25,
             label="mut_wt_rep_vs_mut_wt_conn")

    plt.xlabel("mut_wt_rep vs wt mut conn -log10(p-val)")
    plt.ylabel("Relative frequency")
    plt.legend()

    fig4.savefig("%s_wt_mut_rep_vs_wt_mut_conn_pval_dist.png" % options.output_file_prefix)


    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def buildWT_dict(this_gctx, allele2distil_id, WT_alleles, replicate_null_dist, num_reps):
    """
    {WT_allele:{"wt_rep": med_wt_rep,
                "wt_rep_dist":[]
                "wt_rep_pval": p_val vs null}
    """
    WT_dict = {}
    wt_rep_pvals = []
    wt_allele_ordered = []
    for allele in WT_alleles:
        WT_dict[allele] = {}
        wt_rep_rankpt, rep_rankpts = getSelfConnectivity(this_gctx,
                                                        allele2distil_id[allele],
                                                        num_reps)

        WT_dict[allele]["wt_rep"] = wt_rep_rankpt
        WT_dict[allele]["wt_rep_dist"] = rep_rankpts
        wt_rep_pval = getPairwiseComparisons(rep_rankpts,
                                             replicate_null_dist)
        WT_dict[allele]["wt_rep_pval"] = wt_rep_pval
        wt_rep_pvals.append(wt_rep_pval)
        wt_allele_ordered.append(allele)

    return WT_dict, wt_rep_pvals, wt_allele_ordered

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

#   conn_rankpts = []
#   for i in range(num_reps):
#       for j in range(num_reps):
#           conn_rankpts.append(float(this_gctx.frame[distil_ids1[i]]
#                                                    [distil_ids2[j]]))

#   return numpy.percentile(conn_rankpts, 50), conn_rankpts

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

        # Get replicate null
        # First version. control vs. pert. This approach seemed to not be the
        # best control.
#       rep_null_dist.append(float(this_gctx.frame[random.choice(allele2distil_id[random_control])]
#                                                 [random.choice(allele2distil_id[random_allele])]))

        # Introspect similarity is median of row similarities
        rank_pts = []
        for i in range(1, num_reps):
            rank_pts.append(float(this_gctx.frame[control_distil_ids[0]][control_distil_ids[i]]))

#       for i in range(num_reps):
#           for j in range(i+1,num_reps):
#               rank_pts.append(float(this_gctx.frame[control_distil_ids[i]][control_distil_ids[j]]))

        rep_null_dist.append(numpy.percentile(rank_pts, 50))

#        rep_null_dist.append(float(this_gctx.frame[control_distil_ids[0]]
#                                                  [control_distil_ids[1]]))


        # Get connectivity null
        rank_pts = []
        for i in range(num_reps):
            rank_pts.append(float(this_gctx.frame[allele2distil_id[random_control][0]]
                                                 [allele2distil_id[random_allele][i]]))
#       for i in range(num_reps):
#           for j in range(num_reps):
#               rank_pts.append(float(this_gctx.frame[allele2distil_id[random_control][i]]
#                                                    [allele2distil_id[random_allele][j]]))

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

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

#eVIPpredictfunctions

def get_prediction_1(wt_rep, mut_rep, mut_wt_conn, r_rank, c_rank):
    """
    WT and mut replicate robustness compared using a user-specified rankpoint
    threshold

    Also, mut-wt connectivity is compared when WT and mut are robust
    """
    if wt_rep < r_rank and mut_rep < r_rank:
        return "NI"

    if wt_rep >= r_rank and mut_rep < r_rank:
        return "LOF"

    if wt_rep >= r_rank and mut_rep >= r_rank and mut_wt_conn >= c_rank:
        return "Inert"

    if wt_rep < r_rank and mut_rep >= r_rank:
        return "GOF"

    if wt_rep >= r_rank and mut_rep >= r_rank and mut_wt_conn < c_rank:
        return "GOF"

def get_prediction_2(wt_rep_null_pval, mut_rep_null_pval, mut_wt_conn_null_pval,
                     mut_wt_conn, r_thresh, c_thresh):
    """
    Similar to prediction 1, but using p-value against a null, instead of
    rankpoint thresholds
    """

    if wt_rep_null_pval >= r_thresh and mut_rep_null_pval >= r_thresh:
        return "NI"

    if wt_rep_null_pval < r_thresh and mut_rep_null_pval >= r_thresh:
        return "LOF"

    if wt_rep_null_pval < r_thresh and mut_rep_null_pval < r_thresh and mut_wt_conn_null_pval < c_thresh:
        if mut_wt_conn >= 0:
            return "Inert"
        else:
            return "GOF"

    if wt_rep_null_pval >= r_thresh and mut_rep_null_pval < r_thresh:
        return "GOF"

    if wt_rep_null_pval < r_thresh and mut_rep_null_pval < r_thresh and mut_wt_conn_null_pval >= c_thresh:
        return "GOF"

def get_prediction_3(wt_rep, mut_rep, mut_wt_rep_pval, mut_wt_conn, mut_wt_conn_pval,
                     mut_wt_thresh, mut_wt_rep_diff, c_thresh):
    """
    Instead of checking against a null, a difference in the robustness is used
    """
    if mut_wt_rep_pval < mut_wt_thresh:
        if wt_rep < mut_rep:
            if mut_rep - wt_rep >= mut_wt_rep_diff:
                return "GOF"
        elif wt_rep > mut_rep:
            if wt_rep - mut_rep >= mut_wt_rep_diff:
                return "LOF"

    if mut_wt_conn_pval < c_thresh:
        if mut_wt_conn > 0:
            return "Inert"
        else:
            return "GOF"
    else:
        return "NI"

def get_prediction_4(wt_rep, mut_rep, mut_wt_rep_pval,
                     mut_wt_conn, mut_wt_conn_pval, disting_pval,
                     mut_wt_thresh, mut_wt_rep_diff, c_thresh, disting_thresh):

    if mut_wt_rep_pval < mut_wt_thresh:
        if wt_rep < mut_rep:
            if mut_rep - wt_rep >= mut_wt_rep_diff:
                return "GOF"
        elif wt_rep > mut_rep:
            if wt_rep - mut_rep >= mut_wt_rep_diff:
                return "LOF"

    if mut_wt_conn_pval < c_thresh:
        if mut_wt_conn >= 0:
            if disting_pval < disting_thresh:
                return "COF"
            else:
                return "Inert"
        else:
            return "GOF"
    else:
        return "NI"

def get_prediction_5(wt_rep, mut_rep, mut_wt_rep_pval,
                     mut_wt_conn, mut_wt_conn_pval, disting_pval,
                     mut_wt_thresh, mut_wt_rep_diff, c_thresh, disting_thresh):

    if mut_wt_rep_pval < mut_wt_thresh:
        if wt_rep < mut_rep:
            if mut_rep - wt_rep >= mut_wt_rep_diff:
                return "GOF"
        elif wt_rep > mut_rep:
            if wt_rep - mut_rep >= mut_wt_rep_diff:
                return "LOF"

    if disting_pval < disting_thresh:
        return "COF"

    if mut_wt_conn_pval < c_thresh:
        return "Inert"

    return "NI"

def get_prediction_6(wt_rep, mut_rep, mut_wt_rep_pval,
                     mut_wt_conn, mut_wt_conn_pval, disting_pval,
                     mut_wt_thresh, mut_wt_rep_diff, c_thresh, disting_thresh,
                     conn_null_med):

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
