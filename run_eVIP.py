#!/usr/bin/python
import argparse
import os
import scipy
import eVIP_corr
import eVIP_predict

parser = argparse.ArgumentParser()  # Initiate argument parser object
#from corr
parser.add_argument("-infile", help="Input txt file (filtered and log transformed data")
parser.add_argument("-out_directory", help="path to directory for eVIP output files")
#from compare
parser.add_argument("-sig_info", help = "sig info file with gene information and distil information")
parser.add_argument("-c", help = ".grp file containing allele names of control perturbations. If this file is given, a null will be calculated from these")
parser.add_argument("-r", help = "File explicitly indicating which comparisons to do. Assumes the file has a header and it is ignored. The first column is the reference allele and second column is test allele. If this file is not given, then the reference alleles are assumed to be WT and inferred from the allele names.")
parser.add_argument("-i", help = "Number of iterations to run. DEF=1000")
parser.add_argument("-num_reps", help = "Number of replicates expected for each allele. DEF=3")
parser.add_argument("-ie_col", help = "Name of the column with infection efficiency information. DEF=x_ie_a549")
parser.add_argument("-ie_filter", help = "Threshold for infection efficiency. Any wildtype or mutant alleles having an ie below this threshold, will be removed")
parser.add_argument("-allele_col", help = "Column name that indicates the allele names.DEF=x_mutation_status")
parser.add_argument("-conn_null", help = " Optional file containing connectvity null values from a previous run. Should end in _conn_null.txt")
#from predict
parser.add_argument("-conn_thresh", help = "P-value threshould for connectivity vs null.")
parser.add_argument("-mut_wt_rep_thresh", help = "P-value threshould for comparison of WT and mut robustness")
parser.add_argument("-disting_thresh", help = "P-value threshould that tests if mut and wt reps are indistinguishable from each other")
parser.add_argument("-mut_wt_rep_rank_diff", help = "The minimum difference in median rankpoint WT and mut to consider a difference. DEF=0")
parser.add_argument("-conn_null_med", help = "Median of null connectivity distribution.")
parser.add_argument("-use_c_pval", action ="store_true", help = "Will use corrected p-value instead of raw p-val")
parser.add_argument("-cell_id", help = "Optional: Will only look at signatures from this cell line. Helps to filter sig_info file.")
parser.add_argument("-plate_id", help = "Optional: Will only look at signatures from this plate")

args = parser.parse_args()  # Call command line parser method

#make eVIP output directory
out_dir = args.out_directory
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

#run eVIP_corr.py
print('calculating correlations...')
run_corr = eVIP_corr.api_main(input=args.infile, z_output= out_dir+"/z_scores", sp_output=out_dir+"/spearman_rank_matrix")

# #run eVIP_predict.py
# print('predicting...')
# run_predict = eVIP_predict.api_main(sig_info= args.sig_info,
#         gctx= out_dir+"/spearman_rank_matrix.gct",
#         c= args.c, o= out_dir+"/predict",
#         mut_wt_rep_thresh= args.mut_wt_rep_thresh,
#         disting_thresh=args.disting_thresh, r=args.r,
#         num_reps = args.num_reps , ie_col=args.ie_col,
#         ie_filter=args.ie_filter, i=args.i,
#         conn_thresh=args.conn_thresh,
#         mut_wt_rep_rank_diff= args.mut_wt_rep_rank_diff,
#         conn_null_med = args.conn_null_med, allele_col = args.allele_col,
#         cell_id = args.cell_id, plate_id = args.plate_id,
#         use_c_pval= args.use_c_pval, conn_null = args.conn_null)
#
# print "predict done"
