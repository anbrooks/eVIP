#!/usr/bin/python
import argparse
import os
import urllib2
import errno
import mygene
import scipy
import eVIP_corr
import eVIP_predict
import eVIP_sparkler
import eVIP_viz
from csv import DictWriter

#
# import numpy as np
# import pandas as pd

########
# MAIN #
########
def main(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):

    parser = argparse.ArgumentParser()

    #eVIPP
    parser.add_argument("-eVIPP", action ="store_true", help="Use this option when doing pathway analysis, must also have -file_hsa and -data_file ")
    parser.add_argument("-file_hsa", help="file containing a list of kegg pathway numbers, for running eVIP on several pathways at once ")

    #from corr
    parser.add_argument("--infile", help="Input txt file (filtered and log transformed data)")
    parser.add_argument("-zscore_gct", help="Optional z_score input gct file (use instead of --infile)")
    parser.add_argument("-out_directory",required=True, help="path to directory for eVIP output files")
    #from compare
    parser.add_argument("-sig_info",required=True, help = "sig info file with gene information and distil information")
    parser.add_argument("-c",required=True, help = ".grp file containing allele names of control perturbations. If this file is given, a null will be calculated from these")
    parser.add_argument("-r", required=True, help = "File explicitly indicating which comparisons to do. Assumes the file has a header and it is ignored. The first column is the reference allele and second column is test allele. If this file is not given, then the reference alleles are assumed to be WT and inferred from the allele names.")
    parser.add_argument("-num_reps",required=True, help = "Number of replicates expected for each allele. DEF=3")
    parser.add_argument("-ie_filter", help = "Threshold for infection efficiency. Any wildtype or mutant alleles having an ie below this threshold, will be removed")
    parser.add_argument("-ie_col", help = "Name of the column in the sig_info file with infection efficiency information. DEF=x_ie_a549")
    parser.add_argument("-i", help = "Number of iterations to run. DEF=1000")
    parser.add_argument("-allele_col", help = "Column name in sig_info file that indicates the allele names.DEF=x_mutation_status")
    parser.add_argument("-conn_null", help = " Optional file containing connectvity null values from a previous run. Should end in _conn_null.txt")
    #from predict
    parser.add_argument("-conn_thresh",help = "P-value threshould for connectivity vs null. DEFAULT=0.05")
    parser.add_argument("-mut_wt_rep_thresh", help = "P-value threshould for comparison of WT and mut robustness. DEFAULT=0.05")
    parser.add_argument("-disting_thresh", help = "P-value threshould that tests if mut and wt reps are indistinguishable from each other. DEFAULT=0.05")
    parser.add_argument("-mut_wt_rep_rank_diff", help = "The minimum difference in median rankpoint WT and mut to consider a difference. DEF=0")
    parser.add_argument("-use_c_pval", action ="store_true", help = "Will use corrected p-value instead of raw p-val")
    parser.add_argument("-cell_id", help = "Optional: Will only look at signatures from this cell line. Helps to filter sig_info file.")
    parser.add_argument("-plate_id", help = "Optional: Will only look at signatures from this plate")
    #from sparkler
    parser.add_argument("-ref_allele_mode", action ="store_true", help = "Sparkler+Viz: Instead of organizing plots by gene, will use the wt column to determine what are the reference alleles." )
    parser.add_argument("-x_thresh" , help = "Sparkler: Threshold of significance")
    parser.add_argument("-y_thresh", help = "Sparkler: Threshold of impact direction")
    parser.add_argument("-annotate", action ="store_true", help = "Sparkler: Will add allele labels to points.")
    parser.add_argument("-by_gene_color", help = "Sparkler: File containing labels and colors for gene-centric plot.")
    parser.add_argument("-pdf", help = "Sparkler + Viz: Will print plots in pdf format instead of png.")
    parser.add_argument("-xmin", help = "Sparkler: Min value of x-axis. DEF=0")
    parser.add_argument("-xmax", help = "Sparkler: Max value of x-axis. DEF=4")
    parser.add_argument("-ymin", help = "Sparkler: Min value of y-axis. DEF=-3")
    parser.add_argument("-ymax", help = "Sparkler: Min value of y-axis. DEF=3")
    #from viz
    parser.add_argument("-viz_ymin", help = "Viz: Minimum y-value of rep value. DEF=-100")
    parser.add_argument("-viz_ymax", help = "Viz: Maximum y-value of rep value. DEF=100")
    parser.add_argument("-corr_val", help = "Viz: String used to label the correlation value. DEF= 'row median rankpoints' ")

    global args
    args = parser.parse_args()

    #make eVIP output directory
    global out_dir
    out_dir = args.out_directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    eVIP_dir = out_dir+ "/eVIP_out"
    if not os.path.exists(eVIP_dir):
        os.makedirs(eVIP_dir)

    run_eVIP(args.infile, args.zscore_gct, eVIP_dir, args.sig_info, args.c, args.r, args.num_reps,
                      args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                      args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id, args.ref_allele_mode,
                      args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                      args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)

    #pathway analysis
    if args.eVIPP:
        print "Running eVIPP..."
        hsa_list = open(args.file_hsa, "r")
        global summary_file
        summary_file = open(out_dir + "/eVIPP_summary.txt", "w")
        pathway_list = []
        for line in hsa_list:
          pathway_list.append(line.strip())

        for pathway in pathway_list:
            # make eVIP_out pathway directory
            p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"
            mkdir_p(p_out_directory)

            #running eVIP when the input is zscores
            if args.zscore_gct:

                # make file for the extracted data
                pathway_out_file = open(out_dir + "/" + pathway + "_zscores.gct", "w")
                ensembl_list = kegg_to_ensembl_list(pathway)

                data = open(args.zscore_gct, "r")
                gene_extraction_z(ensembl_list, data, pathway_out_file)

                #zscore path
                zscore_gct = out_dir + "/" + pathway + "_zscores.gct"

                pathway_out_file.close()

                run_eVIP(infile=None, zscore_gct=zscore_gct, out_directory=p_out_directory, sig_info=args.sig_info,
                      c = args.c, r = args.r, num_reps = args.num_reps,
                      ie_filter = args.ie_filter, ie_col = args.ie_col, i = args.i, allele_col = args.allele_col,
                      conn_null = args.conn_null, conn_thresh = args.conn_thresh,
                      mut_wt_rep_rank_diff = args.mut_wt_rep_rank_diff, use_c_pval = args.use_c_pval,
                      cell_id = args.cell_id, plate_id = args.plate_id, ref_allele_mode= args.ref_allele_mode,
                      x_thresh = args.x_thresh, y_thresh = args.y_thresh, annotate = args.annotate,
                      by_gene_color = args.by_gene_color, pdf = args.pdf, xmin = args.xmin,
                      xmax = args.xmax, ymin= args.ymin, ymax = args.ymax, viz_ymin=args.viz_ymin,
                      viz_ymax =args.viz_ymax, corr_val = args.corr_val)


            #running eVIP when the input is expression data
            if args.infile:
                 all_data = open(args.infile, "r")
                 pathway_out_file = open(out_dir + "/" + pathway + "_expression.txt", "w")
                 ensembl_list = kegg_to_ensembl_list(pathway)
                 gene_extraction_infile(ensembl_list, all_data, pathway_out_file)

                 pathway_out_file_txt = out_dir + "/" + pathway + "_expression.txt"
                 pathway_out_file.close()

                 run_eVIP(pathway_out_file_txt, None, p_out_directory, args.sig_info, args.c, args.r, args.num_reps,
                      args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                      args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id, args.ref_allele_mode,
                      args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                      args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)


        summarize(pathway_list)




############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def run_eVIP(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):

    #different sig_gctx for exp an z inputs used in viz

    if args.infile :
        sig_gctx_val = out_directory+ "/z_scores.gct"
    if args.zscore_gct :
        sig_gctx_val = args.zscore_gct




    #run eVIP_corr.py
    print('calculating correlations...')
    run_corr = eVIP_corr.run_main(input=infile,zscore_gct=zscore_gct, out_dir= out_directory)

    #run eVIP_predict.py
    print('predicting...')
    run_predict = eVIP_predict.eVIP_run_main(sig_info= sig_info, o= out_directory+"/predict",c= args.c,r=args.r,
        gctx= out_directory+"/spearman_rank_matrix.gct", conn_thresh=args.conn_thresh, conn_null = args.conn_null,
        allele_col = args.allele_col, ie_filter=args.ie_filter, ie_col=args.ie_col, cell_id = args.cell_id,
        plate_id = args.plate_id, i=args.i, num_reps = args.num_reps, mut_wt_rep_thresh = args.mut_wt_rep_thresh,
        disting_thresh = args.disting_thresh, mut_wt_rep_rank_diff = args.mut_wt_rep_rank_diff, use_c_pval=args.use_c_pval)

    #run eVIP_sparkler.py
    print "making sparkler plots..."
    run_sparkler = eVIP_sparkler.eVIP_run_main(pred_file = out_directory+"/predict.txt", ref_allele_mode=args.ref_allele_mode,
            y_thresh = args.y_thresh , x_thresh = args.x_thresh,
            use_c_pval= args.use_c_pval,annotate=args.annotate, by_gene_color= args.by_gene_color, pdf= args.pdf,
            xmin= args.xmin, xmax = args.xmax, ymin = args.ymin, ymax = args.ymax, out_dir = out_directory+"/sparkler_plots")

    #run eVIP_viz.py
    print "making visualizations..."
    run_viz = eVIP_viz.eVIP_run_main(pred_file= out_directory+"/predict.txt", sig_info = args.sig_info, gctx=out_directory+"/spearman_rank_matrix.gct",
            sig_gctx = sig_gctx_val, ref_allele_mode = args.ref_allele_mode, null_conn = out_directory+"/predict_conn_null.txt",
            out_dir = out_directory+"/viz",ymin = args.viz_ymin, ymax= args.viz_ymax, allele_col = args.allele_col, use_c_pval = args.use_c_pval,
             pdf = args.pdf, cell_id = args.cell_id, plate_id = args.plate_id, corr_val_str= args.corr_val)

    print "eVIP is DONE"


def kegg_to_ensembl_list(hsa_num):
    #this function takes a kegg pathway number, finds the genes present in the pathway,
    #and converts the genes to ensembl ids

    #saving kegg website data
    response = urllib2.urlopen('http://rest.kegg.jp/get/'+hsa_num).read()

    ##parsing the kegg site

    #creating tuple for start and end indices of the genes
    gene_idx_range = tuple()


    #finding the indices

    check_points = ('COMPOUND','REFERENCE')

    for index, item in enumerate(response.split()):
        if item == 'GENE':
            gene_idx_range += (index,)

        if item in check_points:
            gene_idx_range += (index,)
            break

    print "\n"
    print "Pathway:"
    print hsa_num

    #segmenting the data based off the start and end indices
    kegg_items = (response.split()[gene_idx_range[0]:gene_idx_range[1]])

    #getting the gene names
    gene_names = []
    for item in kegg_items:
        if ";" in item:
            item = item.replace(";", "")
            gene_names.append(item)

    print "Genes in pathway:"
    print gene_names

    #converting from official gene symbol to ensembl id
    mg = mygene.MyGeneInfo()
    gene_dict_list = mg.querymany(gene_names, scopes='symbol',fields='ensembl.gene',species='human')

    gene_dicts = []

    for dict in gene_dict_list:
        for key, value in dict.iteritems():
            if key == "ensembl":
                gene_dicts.append(value)

    print "# of gene names:"
    print len(gene_dicts)
    count = 0
    ensembl_list = []
    for d in gene_dicts:
        try:
            value = d.get('gene')
            ensembl_list.append(value)
            count += 1
        except AttributeError:
            for item in d:
                value = item.get('gene')
                ensembl_list.append(value)
                count += 1

    print "# of ensembl ids found:"
    print(count)

    return ensembl_list

def gene_extraction_z(ensembl_list, data, output):
# this function extracts the pathway genes from the RNA-seq data
    matched_ids = []
    header = None
    ncols = None
    ncol_vals = None
    outlines = []

    for line in data:
        line = formatLine(line)
        for id in ensembl_list:
            if line[0].startswith(id):
                matched_ids.append(id)

        ncols = (str(len(line.split())))
        if line.startswith("id"):
            header = line
            ncol_vals = int(ncols)-1

        line = line.split()
        for id in ensembl_list:
            if line[0].startswith(id):
                matched_ids.append(id)
                outline = "\t".join(line) + "\n"
                outlines.append(outline)

    output.write("#1.3\n")
    output.write(str(len(matched_ids)) +"\t" +str(ncol_vals) + "\t" + "0" + "\t" + "0" + "\n")
    output.write(header+"\n")
    for outline in outlines:
        output.write(outline)

    global matched_length_list
    matched_length = str(len(matched_ids))
    ensembl_length = str(len(ensembl_list))
    print matched_length + " of " + ensembl_length + " IDs were found in the data. "

def gene_extraction_infile(ensembl_list, data, output):
    matched_ids = []
    for line in data:
        # line = formatLine(line)
        if line.startswith("#gene_id"):
            output.write(line)
        line = line.split()
        for id in ensembl_list:
            if line[0].startswith(id):
                matched_ids.append(id)
                output.write("\t".join(line) + "\n")

    matched_length = str(len(matched_ids))
    ensembl_length = str(len(ensembl_list))
    print matched_length + " of " + ensembl_length + " IDs were found in the data. "

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def formatLine(line):
    line = line.replace("\r", "")
    line = line.replace("\n", "")
    return line

def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
                yield line.split(delimiter, 1)[1]
            except IndexError:
                continue

def summarize(pathway_list):
    mut_list = []
    pred_list = []

    #getting the list of tested mutations using the first pathway in the list
    with open(out_dir + "/" + pathway_list[0] + "_eVIPP_outputs/predict.txt") as first_pathway:
        for line in first_pathway:
            if line.startswith("gene"):
                continue
            else:
                splitLine = line.split()
                mut_list.append(str(splitLine[1]))

    header = "pathway" + "\t" + ("\t").join(mut_list) + "\t" + "num_genes" + "\n"
    summary_file.write(header)

    for pathway in pathway_list:
        with open(out_dir + "/" + pathway + "_eVIPP_outputs/predict.txt") as this_predict:
            for line in this_predict:
                if line.startswith("gene"):
                    continue
                else:
                    splitLine = line.split()
                    pred_list.append(splitLine[13])

        num_lines_in_zscore = sum(1 for line in open(out_dir + "/" + pathway + "_eVIPP_outputs/z_scores.gct"))

        # summary_file.write(pathway + "\t" + ("\t").join(pred_list)+"\n")

        summary_file.write(pathway + "\t" + ("\t").join(pred_list)+ "\t" + str(num_lines_in_zscore-3)+"\n")


        #clearing predict list
        pred_list = []


#################
# END FUNCTIONS #
#################

if __name__ == "__main__": main()
