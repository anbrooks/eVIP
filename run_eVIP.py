#!/usr/bin/python
import argparse
import os
import errno
import csv
import eVIP_corr
import eVIP_predict
import eVIP_sparkler
import eVIP_viz
import eVIP_compare
import itertools
import rpy2.robjects as robjects
import json
#import eVIPP_sparkler

########
# MAIN #
########
def main(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):

    parser = argparse.ArgumentParser()

    #from corr
    parser.add_argument("--infile", help="Input txt file (filtered and log transformed data).")
    parser.add_argument("-zscore_gct", help="Zscore input gct file (use instead of --infile)")
    parser.add_argument("-out_directory",required=True, help="Path to directory for eVIP output files")
    #from compare
    parser.add_argument("-sig_info",required=True, help = "sig info file with gene information and distil information")
    parser.add_argument("-c",required=True, help = ".grp file containing allele names of control perturbations. If this file is given, a null will be calculated from these")
    parser.add_argument("-r", required=True, help = "File explicitly indicating which comparisons to do. Assumes the file has a header and it is ignored. The first column is the reference allele and second column is test allele. If this file is not given, then the reference alleles are assumed to be WT and inferred from the allele names.")
    parser.add_argument("-num_reps",required=True, help = "Number of replicates expected for each allele. DEF=3")
    parser.add_argument("-ie_filter", help = "Threshold for infection efficiency. Any wildtype or mutant alleles having an ie below this threshold, will be removed")
    parser.add_argument("-ie_col", help = "Name of the column in the sig_info file with infection efficiency information. DEF=x_ie_a549")
    parser.add_argument("-i", help = "Number of iterations to run. DEF=1000")
    parser.add_argument("-allele_col", help = "Column name in sig_info file that indicates the allele names.DEF=x_mutation_status")
    parser.add_argument("-conn_null", help = " Optional file containing connectivity null values from a previous run. Should end in _conn_null.txt")
    #from predict
    parser.add_argument("-conn_thresh",help = "P-value threshold for connectivity vs null. DEFAULT=0.05")
    parser.add_argument("-mut_wt_rep_thresh", help = "P-value threshold for comparison of WT and mut robustness. DEFAULT=0.05")
    parser.add_argument("-disting_thresh", help = "P-value threshold that tests if mut and wt reps are indistinguishable from each other. DEFAULT=0.05")
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
    #eVIPP
    parser.add_argument("-eVIPP", action ="store_true", help="Use this option when doing pathway analysis, must also have JSON file ")
    parser.add_argument("-JSON", help= "JSON file created by create_pathway_JSON.py. Contains dictionary of pathways and the associated ids")
    parser.add_argument("-min_genes", help = "Minimum amount of pathway genes found in data to run eVIPP on. DEF = 5")
    parser.add_argument("-viz_off", action ="store_true", help = "Will not perform eVIP viz step")
    parser.add_argument("-sparkler_off", action ="store_true",help = "Will not perform eVIP sparkler step")

    global args
    args = parser.parse_args()

    #make eVIP output directory
    global out_dir
    out_dir = args.out_directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # running eVIP using JSON (running eVIP using only genes from JSON pathways)
    if args.JSON:
        print "JSON"
        print "Extracting pathway genes from data..."

        #making output for eVIP on the JSON genes
        eVIP_JSON_dir = out_dir + "/eVIP_out_JSON/"
        if not os.path.exists(eVIP_JSON_dir):
            os.makedirs(eVIP_JSON_dir)

        JSON_extracted_loc = out_dir + "/eVIP_out_JSON/JSON_extracted_data.gct"

        #extracting the JSON genes from the data
        JSON_extraction()

        print "running eVIP on pathway genes..."

        #if input is expression data (infile option)
        if args.infile:
            run_eVIP(JSON_extracted_loc, args.zscore_gct, out_dir + "/eVIP_out_JSON/" , args.sig_info, args.c, args.r, args.num_reps,
                 args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                 args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id, args.ref_allele_mode,
                 args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                 args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)

        #if input is z scores (zscore_gct option)
        if args.zscore_gct:
            run_eVIP(args.infile, JSON_extracted_loc, out_dir + "/eVIP_out_JSON/" , args.sig_info, args.c, args.r, args.num_reps,
                 args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                 args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id, args.ref_allele_mode,
                 args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                 args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)


    #running using just the input file (original way to run)
    else:
        eVIP_dir = out_dir + "/eVIP_out"
        if not os.path.exists(eVIP_dir):
            os.makedirs(eVIP_dir)

        run_eVIP(args.infile, args.zscore_gct, eVIP_dir, args.sig_info, args.c, args.r, args.num_reps,
                 args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                 args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id, args.ref_allele_mode,
                 args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                 args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)

    #if using JSON to test pathways
    if args.JSON and args.eVIPP:
        print "Running eVIPP using JSON..."
        summary_file = open(out_dir + "/eVIPP_summary.txt", "w")
        summary_eVIPP_vals = open(out_dir + "/eVIPP_combined_predict_files.txt", "w")

        #getting used pathways (removing pathways in JSON that there isn't data for)
        pway_dict, used_pathways = JSON_pway()

        #running eVIP for each pathway
        JSON_eVIP(used_pathways)

        summarize(used_pathways, summary_file)
        summarize_predict_files(used_pathways, summary_eVIPP_vals)

        #eVIPP sparkler
#
#        print "Making allele pathway sparkler plots..."
#        run_eVIPP_sparkler = eVIPP_sparkler.main(pred_file= out_dir + "/eVIPP_combined_predict_files.txt", ref_allele_mode=args.ref_allele_mode, y_thresh=args.y_thresh, x_thresh=args.x_thresh, use_c_pval=args.use_c_pval,
#                  annotate=args.annotate, by_gene_color=args.by_gene_color, pdf=args.pdf, xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax,
#                  out_dir= out_dir + "/eVIPP_sparkler_plots")


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


    # run eVIP_corr.py
    print('calculating correlations...')
    run_corr = eVIP_corr.run_main(input=infile,zscore_gct=zscore_gct, out_dir= out_directory)

    print('comparing...')
    run_compare = eVIP_compare.run_main(sig_info=sig_info, gctx = out_directory+"/spearman_rank_matrix.gct",
                allele_col = args.allele_col, o= out_directory+"/compare", r = args.r,
             c = args.c, i = args.i, conn_null = args.conn_null, ie_col = args.ie_col,
             ie_filter = args.ie_filter, num_reps = args.num_reps, cell_id = args.cell_id, plate_id = args.plate_id)

    print('predicting...')
    run_predict = eVIP_predict.run_main(i= out_directory+"/compare.txt", o= out_directory+"/predict", conn_thresh=args.conn_thresh,
                mut_wt_rep_thresh=args.mut_wt_rep_thresh, mut_wt_rep_rank_diff=args.mut_wt_rep_rank_diff,
                disting_thresh=args.disting_thresh, use_c_pval=args.use_c_pval)


    if not args.sparkler_off:
        print "making sparkler plots..."
        run_sparkler = eVIP_sparkler.eVIP_run_main(pred_file = out_directory+"/predict.txt", ref_allele_mode=args.ref_allele_mode,
                y_thresh = args.y_thresh , x_thresh = args.x_thresh,
                use_c_pval= args.use_c_pval,annotate=args.annotate, by_gene_color= args.by_gene_color, pdf= args.pdf,
                xmin= args.xmin, xmax = args.xmax, ymin = args.ymin, ymax = args.ymax, out_dir = out_directory+"/sparkler_plots")

    if not args.viz_off:
        print "making visualizations..."
        if args.conn_null:
            null_conn = args.conn_null
        else:
            null_conn = out_directory + "/compare_conn_null.txt"

        run_viz = eVIP_viz.eVIP_run_main(pred_file= out_directory+"/predict.txt", sig_info = args.sig_info, gctx=out_directory+"/spearman_rank_matrix.gct",
                sig_gctx = sig_gctx_val, ref_allele_mode = args.ref_allele_mode, null_conn = null_conn,
                out_dir = out_directory+"/viz",ymin = args.viz_ymin, ymax= args.viz_ymax, allele_col = args.allele_col, use_c_pval = args.use_c_pval,
                 pdf = args.pdf, cell_id = args.cell_id, plate_id = args.plate_id, corr_val_str= args.corr_val)

    print "eVIP is DONE"

def make_adj_compare_file(pathway_list, out_directory,eVIPP_adj_pways_mut_wt_rep_c_pvals_from_compare,eVIPP_adj_pways_mut_wt_conn_null_c_pvals_from_compare,eVIPP_adj_pways_wt_mut_rep_vs_wt_mut_conn_c_pvals_from_compare):

    column_headers = ["gene","mut","mut_rep", "wt_rep", "mut_wt_connectivity",
                      "wt", "cell_line", "mut_wt_rep_pval","mut_wt_conn_null_pval",
                      "wt_mut_rep_vs_wt_mut_conn_pval", "mut_wt_rep_c_pval",
                      "mut_wt_conn_null_c_pval","wt_mut_rep_vs_wt_mut_conn_c_pval"]

    for pathway in pathway_list:
        #new file
        adj_compare = open(out_directory + "/" + pathway + "_eVIPP_outputs/adj_compare.txt", "w")

        #writing header to new file
        file_writer = csv.DictWriter(adj_compare, delimiter="\t", fieldnames=column_headers)
        file_writer.writeheader()

        #opening the original compare file
        compare_file = open(out_directory + "/" +pathway+ "_eVIPP_outputs/compare.txt", "r")
        file_reader = csv.DictReader(compare_file, delimiter="\t")
        n = 0
        m = 0
        l = 0

        for line in file_reader:
            for item in eVIPP_adj_pways_mut_wt_rep_c_pvals_from_compare:
                if pathway in item:
                    line['mut_wt_rep_c_pval'] = item[0][n]
                    n+=1
            for item in eVIPP_adj_pways_mut_wt_conn_null_c_pvals_from_compare:
                if pathway in item:
                    line['mut_wt_conn_null_c_pval'] = item[0][m]
                    m += 1
            for item in eVIPP_adj_pways_wt_mut_rep_vs_wt_mut_conn_c_pvals_from_compare:
                if pathway in item:
                    l += 1
            file_writer.writerow(line)


def run_eVIP_multiple_testing_pt1(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):

    #runs eVIP_corr and compare and returns the p values from compare
    #different sig_gctx for exp an z inputs used in viz

    if args.infile :
        sig_gctx_val = out_directory+ "/z_scores.gct"
    if args.zscore_gct :
        sig_gctx_val = args.zscore_gct

    #run eVIP_corr.py
    print('calculating correlations...')
    run_corr = eVIP_corr.run_main(input=infile,zscore_gct=zscore_gct, out_dir= out_directory)

    #run eVIP_compare.py
    print('comparing...')
    run_compare = eVIP_compare.run_main(sig_info=sig_info, gctx = out_directory+"/spearman_rank_matrix.gct",
                allele_col = args.allele_col, o= out_directory+"/compare", r = args.r,
             c = args.c, i = args.i, conn_null = args.conn_null, ie_col = args.ie_col,
             ie_filter = args.ie_filter, num_reps = args.num_reps, cell_id = args.cell_id, plate_id = args.plate_id)

    #getting the p values from the pathway compare file
    mut_wt_rep_pvals_from_compare = []
    mut_wt_conn_null_pvals_from_compare = []
    wt_mut_rep_vs_wt_mut_conn_pvals_from_compare = []

    num_test_alleles = None

    with open(out_directory+"/compare.txt", "r") as compare_file:
        file_reader = csv.DictReader(compare_file, delimiter="\t")

        for row in file_reader:
            mut_wt_rep_pvals_from_compare.append(row['mut_wt_rep_pval'])
            mut_wt_conn_null_pvals_from_compare.append(row['mut_wt_conn_null_pval'])
            wt_mut_rep_vs_wt_mut_conn_pvals_from_compare.append(row['wt_mut_rep_vs_wt_mut_conn_pval'])


    #counting the number of mutations being tested
    with open(out_directory + "/compare.txt", "r") as compare_file:
        file_reader = csv.DictReader(compare_file, delimiter="\t")
        each_line = list(file_reader)
        num_test_alleles=len(each_line)

    return mut_wt_rep_pvals_from_compare, mut_wt_conn_null_pvals_from_compare, wt_mut_rep_vs_wt_mut_conn_pvals_from_compare,num_test_alleles

def run_eVIP_multiple_testing_pt2(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):

    if args.infile :
        sig_gctx_val = out_directory+ "/z_scores.gct"
    if args.zscore_gct :
        sig_gctx_val = args.zscore_gct

    print('predicting...')
    run_predict = eVIP_predict.run_main(i= out_directory+"/adj_compare.txt", o= out_directory+"/predict", conn_thresh=args.conn_thresh,
                mut_wt_rep_thresh=args.mut_wt_rep_thresh, mut_wt_rep_rank_diff=args.mut_wt_rep_rank_diff,
                disting_thresh=args.disting_thresh, use_c_pval=args.use_c_pval)


    if not args.sparkler_off:
        print "making sparkler plots..."
        run_sparkler = eVIP_sparkler.eVIP_run_main(pred_file = out_directory+"/predict.txt", ref_allele_mode=args.ref_allele_mode,
                y_thresh = args.y_thresh , x_thresh = args.x_thresh,
                use_c_pval= args.use_c_pval,annotate=args.annotate, by_gene_color= args.by_gene_color, pdf= args.pdf,
                xmin= args.xmin, xmax = args.xmax, ymin = args.ymin, ymax = args.ymax, out_dir = out_directory+"/sparkler_plots")

    if not args.viz_off:
        if args.conn_null:
            null_conn = args.conn_null
        else:
            null_conn = out_directory + "/compare_conn_null.txt"

        print "making visualizations..."
        run_viz = eVIP_viz.eVIP_run_main(pred_file= out_directory+"/predict.txt", sig_info = args.sig_info, gctx=out_directory+"/spearman_rank_matrix.gct",
                sig_gctx = sig_gctx_val, ref_allele_mode = args.ref_allele_mode, null_conn = null_conn,
                out_dir = out_directory+"/viz",ymin = args.viz_ymin, ymax= args.viz_ymax, allele_col = args.allele_col, use_c_pval = args.use_c_pval,
                 pdf = args.pdf, cell_id = args.cell_id, plate_id = args.plate_id, corr_val_str= args.corr_val)

    print "Pathway is DONE"


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

def summarize(pathway_list, output_file):
    #creates a simple summary file of just variant predictions for each pathway

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

    header = "pathway" + "\t" + ("\t").join(mut_list) + "\n"
    output_file.write(header)

    for pathway in pathway_list:
        with open(out_dir + "/" + pathway + "_eVIPP_outputs/predict.txt") as this_predict:
            for line in this_predict:
                if line.startswith("gene"):
                    continue
                else:
                    splitLine = line.split()
                    pred_list.append(splitLine[13])

        output_file.write(pathway + "\t" + ("\t").join(pred_list)+"\n")
        #clearing predict list
        pred_list = []

def summarize_predict_files(pathway_list, output_file):
    #combines predict files from each pathway

    # getting header using the first pathway in the list
    with open(out_dir + "/" + pathway_list[0] + "_eVIPP_outputs/predict.txt") as first_pathway:
        for line in first_pathway:
            if line.startswith("gene"):
                output_file.write("Pathway" +"\t"+ line)


    for pathway in pathway_list:
        with open(out_dir + "/" + pathway + "_eVIPP_outputs/predict.txt") as this_predict:
            for line in this_predict:
                if line.startswith("gene"):
                    continue
                else:
                    output_file.write(pathway +"\t"+ line)


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

class ZipExhausted(Exception):
    pass

def izip_longest(*args, **kwds):
    fillvalue = kwds.get('fillvalue')
    counter = [len(args) - 1]
    def sentinel():
        if not counter[0]:
            raise ZipExhausted
        counter[0] -= 1
        yield fillvalue
    fillers = repeat(fillvalue)
    iterators = [chain(it, sentinel(), fillers) for it in args]
    try:
        while iterators:
            yield tuple(map(next, iterators))
    except ZipExhausted:
        pass

def repeat(object, times=None):
    # repeat(10, 3) --> 10 10 10
    if times is None:
        while True:
            yield object
    else:
        for i in xrange(times):
            yield object

def chain(*iterables):
    # chain('ABC', 'DEF') --> A B C D E F
    for it in iterables:
        for element in it:
            yield element

def padj(c_pval_type, num_test_alleles, pathway_list):
    # converting list of lists into one list
    c_pval_type = list(itertools.chain(*c_pval_type))
    #multiple testing corrections
    eVIPP_c_pval_type = robjects.r['p.adjust'](robjects.FloatVector(c_pval_type), "BH")

    new_cpval_list = str(eVIPP_c_pval_type).split()
    fixed_cpval_list = []

    for item in new_cpval_list:
        if item.startswith("["):
            pass
        else:
            fixed_cpval_list.append(item)

    #grouping pvalues from the same pathway together
    grouped_new_cpval_list = list(grouper(fixed_cpval_list, num_test_alleles, fillvalue=None))
    #zipping the pathway with the pvalues
    zipped = zip(grouped_new_cpval_list, pathway_list)

    return zipped

def JSON_extraction():
    all_pway_genes = []

    with open(args.JSON, "rb") as JSON:
        dict = json.load(JSON)
        # creating new dict to remove pathways with None values (no genes were in the pathway)
        new_dict = {}
        new_dict = {k: v for k, v in dict.iteritems() if v is not None}

        # creating list of unique ensembl ids
        for a in new_dict.values():
            for item in a:
                if item in all_pway_genes:
                    continue
                else:
                    all_pway_genes.append(item)

    if args.infile:
        with open(args.infile, "rb") as infile, open(out_dir + "/eVIP_out_JSON/JSON_extracted_data.gct",
                                                     "w") as JSON_extracted:
            matched_ids = []
            for line in infile:
                if line.startswith("#gene_id"):
                    JSON_extracted.write(line)
                sline = line.split()

                for id in all_pway_genes:
                    if sline[0]== id:
                        matched_ids.append(id)
                        JSON_extracted.write("\t".join(sline) + "\n")

    if args.zscore_gct:
        with open(args.zscore_gct, "rb") as zscores, open(out_dir + "/eVIP_out_JSON/JSON_extracted_data.gct", "w") as JSON_extracted:

            matched_ids = []
            outlines = []
            meta_row_names = []
            meta_rows = []

            # saving first line
            first_line = next(zscores)

            # saving second line (#data rows, # data columns, #row metadata fields, #col metadata fields)
            second_line = next(zscores)
            spl = second_line.split()
            ncols = spl[1]
            num_meta_lines= int(spl[3])

            # saving third line (sample names)
            sample_names = next(zscores)

            for n in range(num_meta_lines):
                # getting first element (row name)
                row = (next(zscores).split("\t"))
                meta_rows.append(row)
                meta_row_names.append(row[0])


            for line in zscores:

                sline = line.split()
                for id in all_pway_genes:
                    if sline[0]==id:
                        matched_ids.append(id)
                        outline = "\t".join(sline) + "\n"
                        outlines.append(outline)

            JSON_extracted.write(first_line)
            JSON_extracted.write(str(len(matched_ids)) + "\t" + str(ncols) + "\t" + "0" + "\t" + str(num_meta_lines) + "\n")
            JSON_extracted.write(sample_names)

            for row in meta_rows:
                JSON_extracted.write(str(row[0]) +"\t"+ str("\t".join(row[1:])).strip("\n"))
                JSON_extracted.write("\n")

            for outline in outlines:
                JSON_extracted.write(outline)


def JSON_pway():
    pway_dict = {}
    with open(args.JSON, "rb") as JSON:
        old_dict = json.load(JSON)
        #creating new dict to remove pathways with None values (no genes were in the pathway)
        pway_dict = {k: v for k, v in old_dict.iteritems() if v is not None}

    used_pathways = []

    #for each pathway
    for pathway, genes in pway_dict.iteritems():
        print pathway

        pathway_out_file_txt = out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt"
        p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"
        mkdir_p(p_out_directory)

        pathway_genes = []
        for gene in genes:
            gene = str(gene)
            pathway_genes.append(gene)

        #running eVIP when the input is expression data
        if args.infile:
            with open(out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt", "w+") as pathway_out_file, open(args.infile, "r") as all_data:
                matched_ids = []
                for line in all_data:
                    if line.startswith("#gene_id"):
                        pathway_out_file.write(line)
                    sline = line.split()

                    for id in pathway_genes:
                        if sline[0]==id:
                            matched_ids.append(id)
                            pathway_out_file.write("\t".join(sline) + "\n")

                matched_length = str(len(matched_ids))
                ensembl_length = str(len(pathway_genes))

                print matched_length + " of " + ensembl_length + " IDs were found in the data. "

                percent_matched = (float(len(matched_ids))/(float(len(pathway_genes))))*100
                print "percent matched: " + str(percent_matched)

            if float(matched_length) > (int(args.min_genes) if args.min_genes else 4):
                used_pathways.append(pathway)

            else:
                print "Less than 2 matched genes...skipping eVIP for pathway"

            print "used pathways:"
            print used_pathways

        # running eVIP when the input is expression data
        if args.zscore_gct:
            with open(out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_zscores.gct","w+") as pathway_out_file, open(args.zscore_gct, "r") as all_data:

                matched_ids = []
                outlines = []
                meta_row_names = []
                meta_rows = []

                # saving first line
                first_line = next(all_data)

                # saving second line (#data rows, # data columns, #row metadata fields, #col metadata fields)
                second_line = next(all_data)
                spl = second_line.split()
                ncols = spl[1]
                num_meta_lines = int(spl[3])

                # saving third line (sample names)
                sample_names = next(all_data)

                for n in range(num_meta_lines):
                    # getting first element (row name)
                    row = (next(all_data).split("\t"))
                    meta_rows.append(row)
                    meta_row_names.append(row[0])


                for line in all_data:
                    line = formatLine(line)
                    for id in pathway_genes:
                        if line[0]==id:
                            matched_ids.append(id)

                    sline = line.split()

                    if sline[0] == "id":
                        header = line
                        ncol_vals = int(ncols) - 1

                    for id in pathway_genes:
                        if sline[0]==id:
                            matched_ids.append(id)
                            outline = "\t".join(sline) + "\n"
                            outlines.append(outline)

                pathway_out_file.write(first_line)
                pathway_out_file.write(str(len(matched_ids)) + "\t" + str(ncols) + "\t" + "0" + "\t" + str(num_meta_lines) + "\n")
                pathway_out_file.write(sample_names)

                for row in meta_rows:
                    pathway_out_file.write(str(row[0]) + "\t" + str("\t".join(row[1:])).strip("\n"))
                    pathway_out_file.write("\n")

                for outline in outlines:
                    pathway_out_file.write(outline)

                matched_length = str(len(matched_ids))
                ensembl_length = str(len(pathway_genes))

                print matched_length + " of " + ensembl_length + " IDs were found in the data. "

                percent_matched = (float(len(matched_ids))/(float(len(pathway_genes))))*100
                print "percent matched: " + str(percent_matched)

                if float(matched_length) > (int(args.min_genes) if args.min_genes else 4):
                    used_pathways.append(pathway)

                else:
                    print "Less than 2 matched genes...skipping eVIP for pathway"

    return pway_dict, used_pathways

def JSON_eVIP(used_pathways):

    pways_mut_wt_rep_pvals_from_compare = []
    pways_mut_wt_conn_null_pvals_from_compare = []
    pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare = []

    for pathway in used_pathways:
        pathway_out_file_txt = out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt"
        z_pathway_out_file_txt = out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_zscores.gct"
        p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"
        mkdir_p(p_out_directory)

        if args.infile:
            with open(out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt","r") as pathway_out_file:


                if args.use_c_pval:
                    print "Running with multiple testing..."

                    # running eVIP corr and compare to get pvals
                    mut_wt_rep_pvals_from_compare, mut_wt_conn_null_pvals_from_compare, wt_mut_rep_vs_wt_mut_conn_pvals_from_compare,num_test_alleles \
                        = run_eVIP_multiple_testing_pt1(pathway_out_file_txt, None, p_out_directory, args.sig_info, args.c, args.r, args.num_reps,
                                                        args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                                                        args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id, args.ref_allele_mode,
                                                        args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                                                        args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)

                    # appending values from each pathway to one list
                    pways_mut_wt_rep_pvals_from_compare.append(mut_wt_rep_pvals_from_compare)
                    pways_mut_wt_conn_null_pvals_from_compare.append(mut_wt_conn_null_pvals_from_compare)
                    pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare.append(wt_mut_rep_vs_wt_mut_conn_pvals_from_compare)

                    # if the pathway lists are done being collected
                    if len(pways_mut_wt_rep_pvals_from_compare) == len(used_pathways):
                        eVIPP_mut_wt_rep_pvals = padj(pways_mut_wt_rep_pvals_from_compare, num_test_alleles, used_pathways)

                    if len(pways_mut_wt_conn_null_pvals_from_compare) == len(used_pathways):
                        eVIPP_pways_mut_wt_conn_null_pvals = padj(pways_mut_wt_conn_null_pvals_from_compare,num_test_alleles, used_pathways)

                    if len(pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare) == len(used_pathways):
                        eVIPP_pways_wt_mut_rep_vs_wt_mut_conn_pvals = padj(pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare,num_test_alleles, used_pathways)

                        #adding the new calculated eVIPP pathways to a new compare file
                        make_adj_compare_file(used_pathways, out_dir, eVIPP_mut_wt_rep_pvals,eVIPP_pways_mut_wt_conn_null_pvals,eVIPP_pways_wt_mut_rep_vs_wt_mut_conn_pvals)

                        for pathway in used_pathways:
                            print "pathway running pt2 of eVIP"
                            print pathway

                            pathway_out_file_txt = out_dir + "/" + pathway + "_expression.txt"
                            p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"

                            #getting predictions and finisihin eVIP with the new corrected values
                            run_eVIP_multiple_testing_pt2(pathway_out_file_txt, None, p_out_directory, args.sig_info, args.c,
                                                          args.r, args.num_reps, args.ie_filter, args.ie_col, args.i, args.allele_col,
                                                          args.conn_null, args.conn_thresh,args.mut_wt_rep_rank_diff, args.use_c_pval,
                                                          args.cell_id,args.plate_id, args.ref_allele_mode,args.x_thresh, args.y_thresh,
                                                          args.annotate, args.by_gene_color,args.pdf, args.xmin,args.xmax, args.ymin,
                                                          args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)

                    else:
                        run_eVIP(pathway_out_file_txt, None, p_out_directory, args.sig_info, args.c, args.r,
                                 args.num_reps,
                                 args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                                 args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id,
                                 args.ref_allele_mode,
                                 args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                                 args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)
        if args.zscore_gct:
            with open(out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_zscores.gct","r") as pathway_out_file:
                print pathway

                if args.use_c_pval:
                    print "Running with multiple testing..."

                    # running eVIP corr and compare to get c pvals
                    mut_wt_rep_pvals_from_compare, mut_wt_conn_null_pvals_from_compare, wt_mut_rep_vs_wt_mut_conn_pvals_from_compare, num_test_alleles \
                        = run_eVIP_multiple_testing_pt1(None, z_pathway_out_file_txt, p_out_directory, args.sig_info,
                                                        args.c, args.r, args.num_reps,
                                                        args.ie_filter, args.ie_col, args.i, args.allele_col,
                                                        args.conn_null, args.conn_thresh,
                                                        args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id,
                                                        args.plate_id, args.ref_allele_mode,
                                                        args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color,
                                                        args.pdf, args.xmin,
                                                        args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax,
                                                        args.corr_val)

                    # appending values from each pathway to one list
                    pways_mut_wt_rep_pvals_from_compare.append(mut_wt_rep_pvals_from_compare)
                    pways_mut_wt_conn_null_pvals_from_compare.append(mut_wt_conn_null_pvals_from_compare)
                    pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare.append(
                        wt_mut_rep_vs_wt_mut_conn_pvals_from_compare)

                    # if the pathway lists are done being collected
                    if len(pways_mut_wt_rep_pvals_from_compare) == len(used_pathways):
                        eVIPP_mut_wt_rep_pvals = padj(pways_mut_wt_rep_pvals_from_compare, num_test_alleles, used_pathways)

                    if len(pways_mut_wt_conn_null_pvals_from_compare) == len(used_pathways):
                        eVIPP_pways_mut_wt_conn_null_pvals = padj(pways_mut_wt_conn_null_pvals_from_compare,num_test_alleles, used_pathways)

                    if len(pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare) == len(used_pathways):
                        eVIPP_pways_wt_mut_rep_vs_wt_mut_conn_pvals = padj(pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare,num_test_alleles, used_pathways)

                        #adding the new calculated eVIPP pathways to a new compare file
                        make_adj_compare_file(used_pathways, out_dir, eVIPP_mut_wt_rep_pvals,eVIPP_pways_mut_wt_conn_null_pvals,eVIPP_pways_wt_mut_rep_vs_wt_mut_conn_pvals)

                        for pathway in used_pathways:
                            print "pathway running pt2"
                            print pathway

                            pathway_out_file_txt = out_dir + "/" + pathway + "_expression.txt"
                            p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"

                            # getting predictions and finisihin eVIP with the new corrected values
                            run_eVIP_multiple_testing_pt2(None, z_pathway_out_file_txt, p_out_directory, args.sig_info,
                                                          args.c, args.r, args.num_reps, args.ie_filter, args.ie_col, args.i,
                                                          args.allele_col,args.conn_null, args.conn_thresh, args.mut_wt_rep_rank_diff,
                                                          args.use_c_pval,args.cell_id, args.plate_id, args.ref_allele_mode, args.x_thresh,
                                                          args.y_thresh,args.annotate, args.by_gene_color, args.pdf, args.xmin, args.xmax,
                                                          args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)



                else:
                    run_eVIP(None, z_pathway_out_file_txt, p_out_directory, args.sig_info, args.c, args.r,
                             args.num_reps,
                             args.ie_filter, args.ie_col, args.i, args.allele_col, args.conn_null, args.conn_thresh,
                             args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id,
                             args.ref_allele_mode,
                             args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                             args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)

    return(used_pathways)


#################
# END FUNCTIONS #
#################

if __name__ == "__main__": main()
