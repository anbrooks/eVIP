#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# mutation_impact_viz.py 
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
import csv
import random

from eVIP_compare import getSelfConnectivity, getConnectivity

import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch

import cmap.io.gct as gct

#############
# CONSTANTS #
#############
PRED_TYPE = ["GOF","LOF","COF", "DOM-NEG", "Neutral","NI"]
#DEF_PRED_COL = "Scenario_7__decision_tree_COF"
DEF_PRED_COL = "prediction"

# For jitter plots
WT_RANGE = [9,11]
MUT_RANGE = [19,21]
CONN_RANGE = [29,31]
#NULL_CONN_RANGE = [36,44]

JITTER_XTICKS = [10, 20, 30]
XMAX = 40

DEF_YMIN = -100
DEF_YMAX = 100
DEF_CORR_VAL_STR = "row median rankpoints"

Z_MIN=-10
Z_MAX=10
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
    opt_parser.add_option("--pred_file",
                          dest="pred_file",
                          type="string",
                          help="""File containing the mutation impact
                                  predictions""",
                          default=None)
#   opt_parser.add_option("--col",
#                         dest="pred_col",
#                         type="string",
#                         help="""Prediciton files have predictions based on
#                                 multiple scenarios. The scenario needs to be
#                                 specified because figures will be plotted in
#                                 the order of GOF, LOF, COF,Inert, NI calls. This
#                                 specifies the name of the column that contains
#                                 the prediction. DEF=%s""" % DEF_PRED_COL,
#                         default=DEF_PRED_COL)
    opt_parser.add_option("--sig_info",
                          dest="sig_info",
                          type="string",
                          help="""sig info file with gene information and distil
                                  information""",
                          default=None)
    opt_parser.add_option("--gctx",
                          dest="gctx",
                          type="string",
                          help="GCTX file with correlations",
                          default=None)
    opt_parser.add_option("--sig_gctx",
                          dest="sig_gctx",
                          type="string",
                          help="""GCTX containing signature data. For L1000, this
                                  would  the Z-score data""",
                          default=None)
    opt_parser.add_option("--ref_allele_mode",
                          dest="ref_allele_mode",
                          action="store_true",
                          help="""Instead of organizing plots by gene, will use
                                  the wt column to determine what are the
                                  reference alleles.""",
                          default=False)
    opt_parser.add_option("--null_conn",
                          dest="null_conn",
                          type="string",
                          help="""File of null connectivity values. This file is
                                  given as output from
                                  eVIP_compare.py. The file ends with
                                  conn_null.txt""",
                          default=None)
    opt_parser.add_option("--out_dir",
                          dest="out_dir",
                          type="string",
                          help="Output directory to put figures",
                          default=None)
    opt_parser.add_option("--ymin",
                          dest="ymin",
                          type="int",
                          help="Minimum y-value of rep value. DEF=%d" % DEF_YMIN,
                          default=DEF_YMIN)
    opt_parser.add_option("--ymax",
                          dest="ymax",
                          type="int",
                          help="Maximum y-value of rep value. DEF=%d" % DEF_YMAX,
                          default=DEF_YMAX)
    opt_parser.add_option("--corr_val_str",
                          dest="corr_val_str",
                          type="string",
                          help="String used to label the correlation value. DEF=\"%s\"" % DEF_CORR_VAL_STR,
                          default=DEF_CORR_VAL_STR)
    opt_parser.add_option("--use_c_pval",
                          dest="use_c_pval",
                          action="store_true",
                          help="Use corrected p-val instead of raw pval",
                          default=False)
    opt_parser.add_option("--pdf",
                          dest="pdf",
                          action="store_true",
                          help="Makes figures in pdf format instead of png",
                          default=False)
    opt_parser.add_option("--cell_id",
                          dest="cell_id",
                          type="string",
                          help="""Indicates which cell line. Helps for filtering
                                  sig_info file""",
                          default=None)
    opt_parser.add_option("--plate_id",
                          dest="plate_id",
                          type="string",
                          help="""Indicates which cell line. Helps for filtering
                                  sig_info file""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--pred_file")
#    opt_parser.check_required("--col")
    opt_parser.check_required("--sig_info")
    opt_parser.check_required("--gctx")
    opt_parser.check_required("--null_conn")
    opt_parser.check_required("--out_dir")

    pred_file = open(options.pred_file)
    pred_col = DEF_PRED_COL

    if os.path.exists(options.out_dir):
        out_dir = os.path.abspath(options.out_dir)
    else:
        os.mkdir(options.out_dir)
        out_dir = os.path.abspath(options.out_dir)
        print "Creating output directory: %s" % out_dir 

    pdf = options.pdf
    use_c_pval = options.use_c_pval

    ymin = options.ymin
    ymax = options.ymax

    ref_allele_mode = options.ref_allele_mode
    
    corr_val_str = options.corr_val_str

    cell_id = options.cell_id
    plate_id = options.plate_id

    sig_info = open(options.sig_info)

    null_conn = getNullConnDist(options.null_conn)

#   null_x_vals = []
#   for val in null_conn:
#       null_x_vals.append(random.uniform(NULL_CONN_RANGE[0], NULL_CONN_RANGE[1]))

    this_gctx = gct.GCT(options.gctx)
    this_gctx.read()

    sig_gctx = gct.GCT(options.sig_gctx)
    sig_gctx.read()

    # Process predictions
    # allele2pvals = {allele:[mut vs wt pval, 
    #                         wt vs mut-wt pval,
    #                         mut-wt conn pval]                            
    (gene2wt,
     gene2allele_call,
     gene2num_alleles,
     allele2pvals) = parse_pred_file(pred_file, pred_col, use_c_pval,
                                     ref_allele_mode)
    

    allele2distil_ids = parse_sig_info(sig_info, cell_id, plate_id)

    for gene in gene2wt:

        this_fig = plt.figure()
        this_fig.set_size_inches((gene2num_alleles[gene]+1)*4,
                                  4*3)

        grid_size = (4, gene2num_alleles[gene] + 1)
        
        wt_heatmap_ax = plt.subplot2grid(grid_size, (0,0))
        wt_im = plot_rep_heatmap(wt_heatmap_ax, 
                         this_gctx.frame, 
                         allele2distil_ids[gene2wt[gene]],
                         allele2distil_ids[gene2wt[gene]],
                         gene2wt[gene],
                         ymin, ymax)

        # WT self connectivity
        wt_self, wt_self_row_medians = getSelfConnectivity(this_gctx,
                                                           allele2distil_ids[gene2wt[gene]],        
                                                           len(allele2distil_ids[gene2wt[gene]]))
        
        # Create consistent x values for the wt reps when plotting
        wt_x_vals = []
        for val in wt_self_row_medians:
            wt_x_vals.append(random.randint(WT_RANGE[0], WT_RANGE[1]))

        # Plot color bar on this axis
        plt.colorbar(wt_im, ax=wt_heatmap_ax, shrink=0.7)

        # Plot allele data
        col_counter = 1

        for type in PRED_TYPE:
            for allele in gene2allele_call[gene][type]:

                # CREATE SCATTERPLOT FIGURE
                plot_signatures(pdf, out_dir, 
                                sig_gctx.frame,
                                gene2wt[gene],
                                allele,
                                allele2distil_ids[gene2wt[gene]],
                                allele2distil_ids[allele])

                # PLOT HEATMAP
                this_hm_ax = plt.subplot2grid(grid_size, 
                                             (0, col_counter))
                plot_rep_heatmap(this_hm_ax,
                                 this_gctx.frame,
                                 allele2distil_ids[allele],
                                 allele2distil_ids[allele],
                                 type + " - " + allele,
                                 ymin, ymax)


                # PLOT WT MUT heatmap
                this_wt_mut_ax = plt.subplot2grid(grid_size,
                                                  (1, col_counter))

                plot_rep_heatmap(this_wt_mut_ax,
                                 this_gctx.frame,
                                 allele2distil_ids[gene2wt[gene]],
                                 allele2distil_ids[allele],
                                 gene2wt[gene] + " vs " + allele,
                                 ymin, ymax)

                # PLOT RANKPOINT ROWS
                this_jitter_ax = plt.subplot2grid(grid_size,
                                                  (2, col_counter))

                mut_self, mt_self_row_medians = getSelfConnectivity(this_gctx,
                                                                    allele2distil_ids[allele],
                                                                    len(allele2distil_ids[allele]))
                wt_mut, wt_mut_row_medians = getConnectivity(this_gctx,
                                                             allele2distil_ids[gene2wt[gene]],
                                                             allele2distil_ids[allele],
                                                             len(allele2distil_ids[allele]))                

                plot_jitter(this_jitter_ax,
                            col_counter,
                            wt_x_vals,
                            wt_self_row_medians,
                            mt_self_row_medians,
                            wt_mut_row_medians,
#                            null_x_vals,
#                            null_conn,
                            allele2pvals[allele][0],
                            allele2pvals[allele][1],
                            use_c_pval,
                            ymin, ymax,
                            corr_val_str)
                            

                # Compared to random connectivity
                conn_ax = plt.subplot2grid(grid_size,
                                           (3, col_counter))

                plot_conn(conn_ax,
                          col_counter,
                          null_conn,
                          wt_mut_row_medians,
                          allele2pvals[allele][2],
                          use_c_pval,
                          corr_val_str)

                col_counter += 1
      
        if pdf:  
            this_fig.savefig("%s/%s_impact_pred_plots.pdf" % (out_dir, gene),
                             format="pdf")
        else:
            this_fig.savefig("%s/%s_impact_pred_plots.png" % (out_dir, gene))
        plt.close(this_fig)
			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getNullConnDist(file_name):
    n_file = open(file_name)

    null_conn = []
    for line in n_file:
        line = formatLine(line)
        null_conn.append(float(line))

    return null_conn

def parse_pred_file(pred_file, pred_col, use_c_pval, ref_allele_mode):
    """
    gene2wt,
    gene2allele_call,
    gene2num_alleles
    allele2pvals = {allele:[mut vs wt pval, 
                             wt vs mut-wt pval,
                             mut-wt conn pval]                            
    """
    csv_reader = csv.DictReader(pred_file, delimiter="\t")

    gene2wt = {}
    gene2allele_call = {}
    gene2num_alleles = {}
    allele2pvals = {}
    for row in csv_reader:   
        gene = row["gene"] 
    
        wt = row["wt"]

        if ref_allele_mode:
            gene = wt

        allele = row["mut"]

        if use_c_pval:
            mut_wt_pval = row["mut_wt_rep_c_pval"]
            wt_vs_mut_wt_pval = row["wt_mut_rep_vs_wt_mut_conn_c_pval"]
            mut_wt_conn_pval = row["mut_wt_conn_null_c_pval"]
        else:
            mut_wt_pval = row["mut_wt_rep_pval"]
            wt_vs_mut_wt_pval = row["wt_mut_rep_vs_wt_mut_conn_pval"]
            mut_wt_conn_pval = row["mut_wt_conn_null_pval"]

        pred = row[pred_col]    

        # Set WT allele
        if gene in gene2wt:
            if gene2wt[gene] != wt:
                print "ERROR: Differing WT allele: %s vs %s" % (gene2wt[gene],
                                                                wt)
                sys.exit(1)
        else:
            gene2wt[gene] = wt

        # Add allele
        if gene not in gene2allele_call:
            # Initialize
            gene2allele_call[gene] = {}
            for type in PRED_TYPE:
                gene2allele_call[gene][type] = [] 
            gene2num_alleles[gene] = 0

        gene2allele_call[gene][pred].append(allele)
        gene2num_alleles[gene] += 1

        allele2pvals[allele] = [mut_wt_pval,
                                wt_vs_mut_wt_pval,
                                mut_wt_conn_pval]

    return gene2wt, gene2allele_call, gene2num_alleles, allele2pvals

def parse_sig_info(sig_info, cell_id=None, plate_id=None):
    allele2distil_ids = {}
    
    csv_reader = csv.DictReader(sig_info, delimiter="\t")
    
    for row in csv_reader:   
        if cell_id:
            if cell_id not in row["cell_id"]:
                continue
        if plate_id:
            if plate_id not in row["sig_id"]:
                continue

        allele = row["x_mutation_status"]
        distil_id_field = row["distil_id"]
        distil_ids = distil_id_field.split("|")
        allele2distil_ids[allele] = distil_ids

    return allele2distil_ids
    
def plot_conn(conn_ax, col_counter, null_conn, wt_mut_row_medians,
              conn_pval_text, use_c_pval, corr_val_str):
   
    conn_ax.hist(np.array(null_conn),
                 histtype='stepfilled',
                 normed=True,
                 color='b',
                 alpha=0.25,
                 label="null")
    conn_ax.hist(np.array(wt_mut_row_medians),
                 histtype='stepfilled',
                 normed=True,
                 color='r',
                 alpha=0.25,
                 label="wt_mut_conn")

    conn_ax.set_xlabel(corr_val_str)

    if col_counter == 1:
        conn_ax.set_ylabel("relative frequency")
                
    conn_ax.set_yticklabels([])

    if use_c_pval:
        pval_text = "cP"
    else:
        pval_text = "P"
    conn_ax.text(0.2,0.8,
                "conn_%s:" % pval_text + conn_pval_text,
                 va='bottom',
                 transform = conn_ax.transAxes,
                 size='x-small')
 
def plot_jitter(jitter_ax, col_counter,
                wt_x_vals,
                wt_self_row_medians, 
                mt_self_row_medians,
                wt_mut_row_medians,
#                null_x_vals,
#                null_conn,
                wt_mut_rep_pval_text,
                wt_vs_mut_wt_conn_pval_text,
                use_c_pval,
                ymin, ymax,
                corr_val_str):
    """
    Will mimic a boxplot jitter by plotting y-values with random x-values within
    a short range
    """
    x_vals = []
    y_vals = []
    for i in range(len(wt_self_row_medians)):
        y_vals.append(wt_self_row_medians[i])
        x_vals.append(wt_x_vals[i])
    for val in mt_self_row_medians:
        y_vals.append(val)
        x_vals.append(random.randint(MUT_RANGE[0], MUT_RANGE[1]))
    for val in wt_mut_row_medians:
        y_vals.append(val)
        x_vals.append(random.randint(CONN_RANGE[0], CONN_RANGE[1]))
#   for i in range(len(null_conn)):
#       y_vals.append(null_conn[i])
#       x_vals.append(null_x_vals[i])

    jitter_ax.plot(x_vals, y_vals,
                   'k.',
                   c=((0,0,0,0.25)))

    jitter_ax.set_ylim(ymin,ymax)
    jitter_ax.set_xlim(0, XMAX)
    if col_counter == 1:
        jitter_ax.set_ylabel(corr_val_str)

    jitter_ax.set_xticks(JITTER_XTICKS)
    jitter_ax.set_xticklabels(["wt",
                               "mut",
                               "wt-mut",
                               "random"],
#                               rotation=45,
#                               ha='right',
                               size='x-small')
                                
    # Add p-value text
    if use_c_pval:
        pval_type = "cP"
    else:
        pval_type = "P"
    jitter_ax.text(2,
                   ymin/3,
                   "wt_vs_mut_%s:\n" % pval_type + wt_mut_rep_pval_text,
                   size='x-small')

    jitter_ax.text(XMAX - 2,
                   ymin/3,
                   "wt_vs_conn_%s:\n" % pval_type + wt_vs_mut_wt_conn_pval_text,
                   ha='right',
                   size='x-small')

def plot_rep_heatmap(heatmap_ax, df, distil_ids1, distil_ids2, title, ymin, ymax):
    heatmap_data = df.loc[distil_ids1,distil_ids2]
    dists = distance.squareform(distance.pdist(heatmap_data))
    clusters = sch.linkage(dists, method="average") 
    den = sch.dendrogram(clusters,color_threshold=np.inf, no_plot=True)

    this_im = heatmap_ax.imshow(heatmap_data.ix[den['leaves'],den['leaves']],
                      cmap=plt.cm.bwr,
                      interpolation="nearest",
                      vmin=ymin,
                      vmax=ymax)

    heatmap_ax.get_xaxis().set_visible(False)
    heatmap_ax.get_yaxis().set_visible(False)
    
    heatmap_ax.set_title(title, size='x-small')

    return this_im

def plot_signatures(pdf, out_dir, sig_gctx_frame, wt_allele, mut_allele,
                    wt_distil_ids, mut_distil_ids):
    num_reps = len(wt_distil_ids)
    this_fig = plt.figure()
    this_fig.set_size_inches(((4*num_reps)*2), (4*num_reps)*2)

    grid_size = (num_reps*2, num_reps*2)

    all_distil_ids = wt_distil_ids + mut_distil_ids
   
    for i in range(num_reps * 2):
        for j in range(i,num_reps*2):
            this_ax = plt.subplot2grid(grid_size, (i,j))
            if i == j:
                if i < num_reps:
                    this_ax.text(0.25,0.5,
                                 wt_allele,
                                 size='large') 
                else:
                    this_ax.text(0.25,0.5,
                                 mut_allele,
                                 size='large',
                                 color='red') 
                continue

            # linear fit to data
            fit = np.polyfit(sig_gctx_frame.loc[:,all_distil_ids[i]],
                             sig_gctx_frame.loc[:,all_distil_ids[j]],
                             deg=1)

            this_ax.plot(sig_gctx_frame.loc[:,all_distil_ids[i]],
                         sig_gctx_frame.loc[:,all_distil_ids[j]],
                         'k.',
                         c=(0,0,0,0.1))

            x_vals = np.arange(Z_MIN, Z_MAX)

            # linear fit plot
            this_ax.plot(x_vals,
                         fit[0]*x_vals + fit[1],
                         "-",c=(1,0,0,0.5))

            # x=y plot
            this_ax.plot(x_vals, x_vals,
                         "--",c=(0,0,0,0.5))

            this_ax.set_xlim(Z_MIN,Z_MAX)
            this_ax.set_ylim(Z_MIN,Z_MAX)
    
    if pdf:
        this_fig.savefig("%s/%s_%s_scatter_plots.pdf" % (out_dir, wt_allele, mut_allele),
                         format="pdf")
    else:
        this_fig.savefig("%s/%s_%s_scatter_plots.png" % (out_dir, wt_allele, mut_allele))

    plt.close(this_fig)
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
