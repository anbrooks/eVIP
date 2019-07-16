#!/usr/bin/python

import sys
import optparse
import os
import pdb
import csv
import random
import math

from eVIP_compare import getSelfConnectivity, getConnectivity
from eVIP_predict import max_diff

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.markers as mmarkers
from matplotlib.colors import colorConverter
import numpy as np
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch

import cmap.io.gct as gct

#############
# CONSTANTS #
#############
DIFF_THRESH = 5
INFINITY = 10
DEF_PRED_COL = "prediction"

#RGB CODES
#GOF_COL = "#e31a1c"
GOF_COL = colorConverter.to_rgba("#ca0020", 1)
COF_PLUS_COL = colorConverter.to_rgba("#c2a5cf", 1)
COF_COL = colorConverter.to_rgba("#6a3d9a", 1)
COF_MINUS_COL = colorConverter.to_rgba("#c2a5cf", 1)
#LOF_COL = "#1f78b4"
LOF_COL = colorConverter.to_rgba("#0571b0", 1)
EQ_COL = colorConverter.to_rgba("#5aae61", 1)
INERT_COL = colorConverter.to_rgba("#000000", 1)
NI_COL = colorConverter.to_rgba("#ffffff", 1)

MAIN_MARKER = 'o'
NEG_MARKER = 'o'


XMIN=0
XMAX=4
YMIN=-3
YMAX=3

MARKER_SIZE= 300
SPARKLER_LINEWIDTH = 3
ALL_SPARKLER_LINEWIDTH = 2

THRESH_LS = "dotted"
#################
# END CONSTANTS #
#################

########
# MAIN #
########

def main(pred_file=None, ref_allele_mode=None, y_thresh=None, x_thresh=None, use_c_pval=None,
                  annotate=None, by_gene_color=None, pdf=None, xmin=None, xmax=None, ymin=None, ymax=None,
                  out_dir=None):

    x_thresh = float(x_thresh)
    y_thresh = float(y_thresh)

    # setting default values
    xmin = float(xmin) if xmin != None else float(0.0)
    xmax = float(xmax) if xmax != None else float(4.0)
    ymin = float(ymin) if ymin != None else float(-3.0)
    ymax = float(ymax) if ymax != None else float(3.0)

    if pdf:
        format = "pdf"
    else:
        format = "png"

    if os.path.exists(out_dir):
        out_dir = os.path.abspath(out_dir)
    else:
        os.mkdir(out_dir)
        out_dir = os.path.abspath(out_dir)
        print "Creating output directory: %s" % out_dir

    x_thresh = getNegLog10(x_thresh, xmax)
    y_thresh = getNegLog10(y_thresh, ymax)
    pred_col = DEF_PRED_COL


    allele2type = None
    if by_gene_color:
        allele2type = parseGeneColor(by_gene_color)


    (allele2mut_wt,
    allele2mut_wt_rep_p,
    allele2neg_log_p,
    allele2diff_score,
    allele2gene,
    allele2col,
    allele2markerstyle) = allele_parse_eVIPP_combined_pred_file(pred_file, x_thresh, y_thresh,
                                         pred_col, use_c_pval, allele2type, ref_allele_mode, xmax, ymax)

    all_mut_wt = []
    all_mut_wt_rep_p = []
    all_neg_log_p = []
    all_col = []
    all_diff_score = []
    all_markerstyle = []
    recs = []
    predictions = ["GOF", "LOF", "COF", "Neutral"]


    for allele, value in allele2gene.iteritems():

        this_fig = plt.figure()
        ax = this_fig.add_subplot(111)

        all_mut_wt.extend(allele2mut_wt[allele])
        all_mut_wt_rep_p.extend(allele2mut_wt_rep_p[allele])
        all_neg_log_p.extend(allele2neg_log_p[allele])
        all_col.extend(allele2col[allele])
        all_diff_score.extend(allele2diff_score[allele])
        all_markerstyle.extend(allele2markerstyle[allele])


        if allele not in allele2type:
            allele2type[allele] = "UNKN"


            (main_markers,
             neg_markers) = split_data(allele2markerstyle[allele],
                                        allele2neg_log_p[allele],
                                       allele2mut_wt_rep_p[allele],
                                       allele2col[allele])

            try:
                plt.scatter(main_markers["x"],
                            main_markers["y"],
                            s=MARKER_SIZE,
                            c=main_markers["col"],
                            marker=MAIN_MARKER,
                            edgecolors="none",
                            linewidth=0)
            except:
                pdb.set_trace()

            plt.scatter(neg_markers["x"],
                        neg_markers["y"],
                        s=MARKER_SIZE,
                        c=neg_markers["col"],
                        marker=NEG_MARKER,
                        linewidth=4)

            for i in range(len(allele2neg_log_p[allele])):
                this_col = allele2col[allele][i]
                if this_col == colorConverter.to_rgba("#ffffff", 1):  # white
                    this_col = "black"
                plt.plot([0, allele2neg_log_p[allele][i]],
                         [0, allele2mut_wt_rep_p[allele][i]],
                         color=this_col,
                         linewidth=SPARKLER_LINEWIDTH)

            if annotate:
                for i in range(len(allele2gene[allele])):
                    ax.annotate(allele2gene[allele][i],
                                (allele2neg_log_p[allele][i],
                                 allele2mut_wt_rep_p[allele][i]),
                                textcoords='data')

            plt.axvline(x= x_thresh, color="grey", ls=THRESH_LS)

            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

            if use_c_pval:
                ax.set_xlabel("-log10(corrected p-val)")
            else:
                ax.set_xlabel("-log10(p-val)")
            ax.set_ylabel("impact direction score")

            ax.text(100, -7, "MUT robustness < WT robustness", fontsize="small",
                    ha="center")
            ax.text(100, 7, "MUT robustness > WT robustness", fontsize="small",
                    ha="center")

            predictions = ["GOF", "LOF", "COF", "Neutral"]
            colors = [GOF_COL, LOF_COL, COF_COL, "black"]

            for i in range(len(colors)):
                recs.append(mpatches.Rectangle((0, 0), 1, 1, fc=colors[i]))

            this_fig.savefig("%s/%s_spark_plots.%s" % (out_dir, allele, format),
                             format=format)

            plt.close(this_fig)


        # final plot
        this_fig = plt.figure()
        ax = this_fig.add_subplot(111)

        # Add lines
        for i in range(len(all_diff_score)):
            this_col = all_col[i]
            if this_col == "white":
                this_col = "black"
            plt.plot([0, all_neg_log_p[i]],
                     [0, all_mut_wt_rep_p[i]],
                     color=this_col,
                     alpha=0.25,
                     linewidth=ALL_SPARKLER_LINEWIDTH)

        (main_markers,
         neg_markers) = split_data(all_markerstyle,
                                   all_neg_log_p,
                                   all_mut_wt_rep_p,
                                   all_col)

        plt.scatter(main_markers["x"],
                    main_markers["y"],
                    s=MARKER_SIZE,
                    c=main_markers["col"],
                    marker=MAIN_MARKER,
                    edgecolors="none",
                    linewidth=0)

        plt.scatter(neg_markers["x"],
                    neg_markers["y"],
                    s=MARKER_SIZE,
                    c=neg_markers["col"],
                    marker=NEG_MARKER,
                    linewidth=4)

        plt.axvline(x=x_thresh, color="grey", ls=THRESH_LS)

        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

        if use_c_pval:
            ax.set_xlabel("-log10(corrected p-val)")
        else:
            ax.set_xlabel("-log10(p-val)")
        ax.set_ylabel("impact direction score")

        this_fig.savefig("%s/all_spark_plots.%s" % (out_dir, format), format=format)

        plt.close(this_fig)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def getNegLog10(p_val, max_val):
    if p_val == 0.0:
        return max_val

    return -math.log(p_val, 10)

def formatAllele(allele):
    allele_elems = allele.split(".")[1:]

    return ".".join(allele_elems)

def parseGeneColor(by_gene_color_file_name):
    """
    Returns
    gene2label
    """
    gene_file = open(by_gene_color_file_name)

    gene2label = {}

    csv_reader = csv.DictReader(gene_file, delimiter="\t")
    for row in csv_reader:
        gene2label[row["gene"]] = row["label"]

    return gene2label


def allele_parse_eVIPP_combined_pred_file(pred_file, x_thresh, y_thresh, pred_col, use_c_pval, allele2type, ref_allele_mode, xmax, ymax):

    #also making dict for alleles so plots can be done per allele not per gene

    allele2mut_wt = {}
    allele2mut_wt_rep_p = {}
    allele2neg_log_p = {}
    allele2diff_score = {}
    allele2gene= {}
    allele2col = {}
    allele2markerstyle = {}

    with open(pred_file) as csvfile:

        csv_reader = csv.DictReader(csvfile, delimiter="\t")
        l = list(csv_reader)

        for row in l:

            gene = row["gene"]
            wt_ref = row["wt"]
            allele = row["mut"]
            mut_rep = float(row["mut_rep"])
            wt_rep = float(row["wt_rep"])
            pred = row[pred_col]

            # Skip NI
            if pred == "NI":
                continue

            if ref_allele_mode:
                gene = wt_ref

            mut_wt_conn = float(row["mut_wt_connectivity"])

            if use_c_pval:
                mut_wt_rep_pval = getNegLog10(float(row["mut_wt_rep_c_pval"]), ymax)

                mut_wt_conn_pval = getNegLog10(float(row["mut_wt_conn_null_c_pval"]), xmax)
                impact_pval = getNegLog10(float(row["wt_mut_rep_vs_wt_mut_conn_c_pval"]), xmax)
            else:

                mut_wt_rep_pval = getNegLog10(float(row["mut_wt_rep_pval"]), ymax)
                mut_wt_conn_pval = getNegLog10(float(row["mut_wt_conn_null_pval"]), xmax)
                impact_pval = getNegLog10(float(row["wt_mut_rep_vs_wt_mut_conn_pval"]), xmax)


            # Initiatite dictionaries

            if allele not in allele2mut_wt:
                #initiate dictionaries
                allele2mut_wt[allele] = []
                allele2mut_wt_rep_p[allele] = []
                allele2neg_log_p[allele] = []
                allele2diff_score[allele] = []
                allele2gene[allele] = []
                allele2col[allele] = []
                allele2markerstyle[allele] = []

            #add gene
            allele2gene[allele].append(gene)

            # MUT - WT
            allele2mut_wt[allele].append(mut_rep - wt_rep)
            allele2neg_log_p[allele].append(impact_pval)

            if mut_rep >= wt_rep:
                allele2mut_wt_rep_p[allele].append(mut_wt_rep_pval)
            else:
                allele2mut_wt_rep_p[allele].append(-mut_wt_rep_pval)

            diff = max_diff(wt_rep, mut_rep, mut_wt_conn)
            allele2diff_score[allele].append(diff)

            # Other features based on significance
            if pred == "GOF":
                allele2col[allele].append(GOF_COL)
                allele2markerstyle[allele].append(MAIN_MARKER)
            elif pred == "LOF":
                allele2col[allele].append(LOF_COL)
                allele2markerstyle[allele].append(MAIN_MARKER)
            elif pred == "Neutral":
                allele2col[allele].append(INERT_COL)
                allele2markerstyle[allele].append(MAIN_MARKER)
            elif pred == "NI":
                allele2col[allele].append(NI_COL)
                allele2markerstyle[allele].append(MAIN_MARKER)
            elif pred == "COF":
                allele2col[allele].append(COF_COL)
                allele2markerstyle[allele].append(MAIN_MARKER)
            elif pred == "DOM-NEG":
                if mut_wt_rep_pval < y_thresh:
                    allele2col[allele].append(COF_COL)
                else:
                    if mut_rep >= wt_rep:
                        allele2col[allele].append(GOF_COL)
                    else:
                        allele2col[allele].append(LOF_COL)
                allele2markerstyle[allele].append(NEG_MARKER)

            else:
                allele2col[allele].append("yellow")
                allele2markerstyle[allele].append(MAIN_MARKER)

    return allele2mut_wt, allele2mut_wt_rep_p, allele2neg_log_p, allele2diff_score, allele2gene, allele2col, allele2markerstyle


def split_data(marker_list, x_list, y_list, col_list):
    """
    Returns dictioaries for x, y, and colors
    """
    main_markers = {"x":[],
                    "y":[],
                    "col":[]}

    neg_markers = {"x":[],
                   "y":[],
                   "col":[]}

    for i in range(len(marker_list)):
        if marker_list[i] == MAIN_MARKER:
            main_markers["x"].append(x_list[i])
            main_markers["y"].append(y_list[i])
            main_markers["col"].append(col_list[i])
        else:
            neg_markers["x"].append(x_list[i])
            neg_markers["y"].append(y_list[i])
            neg_markers["col"].append(col_list[i])

    return main_markers, neg_markers
#################
# END FUNCTIONS #
################
if __name__ == "__main__": main()
