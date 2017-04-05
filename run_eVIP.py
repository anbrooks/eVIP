#!/usr/bin/python
import argparse
import os
import scipy
import eVIP_corr
import eVIP_predict


########
# MAIN #
########
def main():

    parser = argparse.ArgumentParser()  # Initiate argument parser object
    #from corr
    parser.add_argument("--infile",required=True, help="Input txt file (filtered and log transformed data")
    parser.add_argument("-out_directory",required=True, help="path to directory for eVIP output files")
    #from compare
    parser.add_argument("-sig_info",required=True, help = "sig info file with gene information and distil information")
    parser.add_argument("-c",required=True, help = ".grp file containing allele names of control perturbations. If this file is given, a null will be calculated from these")
    parser.add_argument("-r", required=True, help = "File explicitly indicating which comparisons to do. Assumes the file has a header and it is ignored. The first column is the reference allele and second column is test allele. If this file is not given, then the reference alleles are assumed to be WT and inferred from the allele names.")
    parser.add_argument("-i", help = "Number of iterations to run. DEF=1000")
    parser.add_argument("-num_reps",required=True, help = "Number of replicates expected for each allele. DEF=3")
    parser.add_argument("-ie_col", help = "Name of the column with infection efficiency information. DEF=x_ie_a549")
    parser.add_argument("-ie_filter",required=True, help = "Threshold for infection efficiency. Any wildtype or mutant alleles having an ie below this threshold, will be removed")
    parser.add_argument("-allele_col", help = "Column name that indicates the allele names.DEF=x_mutation_status")
    parser.add_argument("-conn_null", help = " Optional file containing connectvity null values from a previous run. Should end in _conn_null.txt")
    #from predict
    parser.add_argument("-conn_thresh",help = "P-value threshould for connectivity vs null. DEFAULT=0.05")
    parser.add_argument("-mut_wt_rep_thresh", help = "P-value threshould for comparison of WT and mut robustness. DEFAULT=0.05")
    parser.add_argument("-disting_thresh", help = "P-value threshould that tests if mut and wt reps are indistinguishable from each other. DEFAULT=0.05")
    parser.add_argument("-mut_wt_rep_rank_diff", help = "The minimum difference in median rankpoint WT and mut to consider a difference. DEF=0")
    parser.add_argument("-conn_null_med", help = "Median of null connectivity distribution.")
    parser.add_argument("-use_c_pval", action ="store_true", help = "Will use corrected p-value instead of raw p-val")
    parser.add_argument("-cell_id", help = "Optional: Will only look at signatures from this cell line. Helps to filter sig_info file.")
    parser.add_argument("-plate_id", help = "Optional: Will only look at signatures from this plate")
    #from sparkler





    args = parser.parse_args()  # Call command line parser method

    #make eVIP output directory
    out_dir = args.out_directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #run eVIP_corr.py
    print('calculating correlations...')
    run_corr = eVIP_corr.api_main(input=args.infile, z_output= out_dir+"/z_scores", sp_output=out_dir+"/spearman_rank_matrix")

    #run eVIP_predict.py
    print('predicting...')
    run_predict = eVIP_predict.eVIP_run_main(sig_info= args.sig_info,
        o= out_dir+"/predict",c= args.c,r=args.r,
        gctx= out_dir+"/spearman_rank_matrix.gct",
        conn_thresh=args.conn_thresh,
        conn_null = args.conn_null,
        allele_col = args.allele_col,
        ie_filter=args.ie_filter,
        ie_col=args.ie_col,
        cell_id = args.cell_id,
        plate_id = args.plate_id,
        i=args.i,
        num_reps = args.num_reps,
        mut_wt_rep_thresh = args.mut_wt_rep_thresh,
        disting_thresh = args.disting_thresh,
        mut_wt_rep_rank_diff = args.mut_wt_rep_rank_diff,
        conn_null_med= args.conn_null_med)

    print "predictions done"

    #run sparkler
    print "making sparkler plots..."
    #run_sparkler = eVIP_sparkler.main()

    #run viz
    print "making visualizations..."
    #run_viz = eVIP_viz.main()

    print "DONE"

############
# END_MAIN #
############

if __name__ == "__main__": main()
