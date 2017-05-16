#!/usr/bin/python
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import rpy2.robjects as robjects
import os

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-input", help="Input txt file (filtered and log transformed data")
    parser.add_argument("-zscore_gct", help="Z-score input gct file (spearman rank correlations will be calculated from this file)")
    parser.add_argument("-out_dir", help = "Directory for output files")

    args = parser.parse_args()

    #if the input is gene expression data
    if args.input != None:
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)

        input_file = open(args.input, "r")
        z_output_file = open(args.out_dir+"z_scores.gct", "w+")
        z_output_file.writelines("#1.3" + "\n")
        sp_output = open(args.out_dir + "spearman_rank_matrix.gct", "w")
        sp_output.writelines("#1.3" + "\n")

        #formatting header and gct files
        for line in open(args.input,"r"):
            line = formatLine(line)
            ncols = (str(len(line.split())))
            if line.startswith("#"):
                header = line
                header += "\n"
                ncol_vals = int(ncols) - 1
                line_count = getLineCount(input_file)
                z_output_file.write(str(line_count-1) + "\t" + str(ncol_vals) + "\t" + "0" + "\t" + "0" + "\n")
                sp_output.write(str(ncol_vals) + "\t" + str(ncol_vals) + "\t" + "0" + "\t" + "0" + "\n")
                header = header.replace("#gene_id", "id")
                sp_output.write(header)
                z_output_file.write(header)
                continue

        #calculating z_scores
        z_outlines = calcZscore(open(args.input, "r"))

        #writing z_scores to output file
        for outline in z_outlines:
            z_output_file.write(outline)

        #importing zscores as matrix
        file_path = args.out_dir+"z_scores.gct"
        data = pd.read_csv(file_path, delimiter="\t", skiprows=2, index_col=False)
        column_name = 'id'
        data = data.drop(column_name, axis=1)

        #calculating spearman rank correlation
        sp_matrix = (stats.spearmanr(data))

        #writing matrix to output file
        matrixToFile(sp_matrix[0], sp_output, header)

    #if the input is zscores
    if args.zscore_gct:
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)

        zscore_input = open(args.zscore_gct, "r")

        #importing zscores as a matrix
        data = np.loadtxt(strip_first_col(args.zscore_gct, "\t"), skiprows=3)

        #calculating spearman rank correlation
        sp_matrix = (stats.spearmanr(data))

        #beggining gct file format
        sp_output = open(args.out_dir+"/spearman_rank_matrix.gct", "w")
        sp_output.writelines("#1.3" + "\n")

        header = None
        ncols = None

        #getting the header and creating gct format
        for line in zscore_input:
            line = formatLine(line)
            ncols = (str(len(line.split())))
            if line.startswith("id"):
                header = line
                header += "\n"
                ncol_vals = int(ncols) - 1
                sp_output.write(str(ncol_vals) + "\t" + str(ncol_vals) + "\t" + "0" + "\t" + "0" + "\n")
                sp_output.write(header)


        # writing matrix to output file
        matrixToFile(sp_matrix[0], sp_output, header)

def run_main(input=None, zscore_gct=None, out_dir=None):
    # if the input is gene expression data
    if input != None:

        with open(input,"r") as input_file:
            z_output_file = open(out_dir + "/z_scores.gct", "w+")
            sp_output = open(out_dir + "/spearman_rank_matrix.gct", "w")
            sp_output.writelines("#1.3" + "\n")
            z_output_file.writelines("#1.3" + "\n")
            for line in input_file:
                line = formatLine(line)
                ncols = (str(len(line.split())))
                if line.startswith("#gene_id"):
                    header = line
                    header += "\n"
                    ncol_vals = int(ncols) - 1
                    z_output_file.write(str(getLineCount(input_file)) + "\t" + str(ncol_vals) + "\t" + "0" + "\t" + "0" + "\n")
                    sp_output.write(str(ncol_vals) + "\t" + str(ncol_vals) + "\t" + "0" + "\t" + "0" + "\n")
                    header = header.replace("#gene_id", "id")
                    sp_output.write(header)
                    z_output_file.write(header)

            z_outlines = calcZscore(open(input, "r"))

            for outline in z_outlines:
                z_output_file.write(outline)

        file_path = out_dir + "/z_scores.gct"

        z_output_file.close()

        data = pd.read_csv(file_path, delimiter="\t", skiprows=2, index_col=False)

        column_name = 'id'
        data = data.drop(column_name, axis=1)
        sp_matrix = (stats.spearmanr(data))

        #writing matrix to output file
        matrixToFile(sp_matrix[0], sp_output, header)

        sp_output.close()


    # if the input is zscores
    elif zscore_gct != None:

        # zscore_copy = out_dir+"/z_scores.gct"
        # shutil.copy(zscore_gct,zscore_copy)

        sp_output = open(out_dir + "/spearman_rank_matrix.gct", "w")
        sp_output.writelines("#1.3" + "\n")

        # data = pd.read_csv(zscore_copy, delimiter="\t", skiprows=2, index_col=False)
        # column_name = 'id'
        # data = data.drop(column_name, axis=1)
        # sp_matrix = (stats.spearmanr(data))

        zscores = open(zscore_gct, "r")
        ncols = None
        header = None

        #getting the header and creating gct format
        for line in zscores:
            line = formatLine(line)
            ncols = (str(len(line.split())))
            if line.startswith("id") :
                header = line
                header += "\n"
                ncol_vals = int(ncols) - 1
                sp_output.write(str(ncol_vals) + "\t" + str(ncol_vals) + "\t" + "0" + "\t" + "0" + "\n")
                sp_output.write(header)
                continue

        zscores.close()

        data = pd.read_csv(zscore_gct, delimiter="\t", skiprows=2, index_col=False)
        column_name = 'id'
        data = data.drop(column_name, axis=1)

        sp_matrix = (stats.spearmanr(data))


        # writing matrix to output file
        matrixToFile(sp_matrix[0], sp_output, header)
        sp_output.close()


#############
# FUNCTIONS #
#############
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getLineCount(file):
    count = 0
    for line in file.xreadlines(): count += 1
    return count

def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
                yield line.split(delimiter, 1)[1]
            except IndexError:
                continue



def matrixToFile(matrix, output_file, header):
    n=1
    for line in matrix:
        sp_line = np.ndarray.tolist(line)
        header_list = header.split("\t")
        sp_line.insert(0, header_list[n])
        for i, value in enumerate(sp_line):
            output_file.write(str(value).replace("\n", ""))
            if i != len(sp_line)-1 :
                output_file.write("\t")
            else:
                output_file.write("\n")
        n += 1

def calcZscore(input_file):
    z_outline = []
    z_outlines = []
    for line in input_file:
        if line.startswith("#"):
            continue

        line = formatLine(line)
        ncols = (str(len(line.split())))
        lineList = line.split("\t")

        vals = []
        in_vals = []
        for item in lineList[1:]:
            val = float(item)
            in_vals.append(val)
            vals.append(val)

        median = robjects.r['mean'](robjects.FloatVector(vals))[0]
        mad = robjects.r['sd'](robjects.FloatVector(vals))[0]

        z_vals = []
        t_vals = []

        for v in in_vals:
            if mad == 0:
                z_val = (v - median) / CLOSE_TO_ZERO
            else:
                z_val = (v - median) / mad

            z_vals.append(z_val)
            t_vals.append(z_val)

        z_vals_str = []
        for v in t_vals:
            z_vals_str.append("%.4f" % v)

        z_outline = "\t".join(lineList[:1]) + "\t"
        z_outline += "\t".join(z_vals_str)
        z_outline += "\n"


        z_outlines.append(z_outline)

    return z_outlines


#################
# END FUNCTIONS #
#################
if __name__ == "__main__": main()

