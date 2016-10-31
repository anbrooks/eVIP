import sys
import optparse
import numpy as np
from scipy import stats
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

#############
# CONSTANTS #
#############

CLOSE_TO_ZERO = 0.0001

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
    opt_parser.add_option("-i",
                          dest="input",
                          type="string",
                          help="""input txt file with filtered gene expression data""",
                          default=None)
    opt_parser.add_option("--z_out",
                          dest="z_output",
                          type="string",
                          help="""name of output file containing z_scores""",
                          default=None)
    opt_parser.add_option("--sp_out",
                          dest="sp_output",
                          type="string",
                          help="""name of output file containing spearman correlation""",
                          default=None)

    (options, args) = opt_parser.parse_args()
    opt_parser.check_required("-i")
    opt_parser.check_required("--z_out")
    opt_parser.check_required("--sp_out")

    input_file = open(options.input, "r")
    z_output= open(options.z_output+".txt", "w+")
    sp_output = open(options.sp_output+".txt", "w")

    t_vals=[]
    z_outline = []
    z_outlines = []

    sp_output.writelines("#1.3"+"\n")




    header = None
    for line in input_file:
        line = formatLine(line)
        ncols = (str(len(line.split())))
        if line.startswith("#"):
            header = line
            header += "\n"
            ncol_vals = int(ncols)-1
            print ncol_vals
            sp_output.write(str(ncol_vals)+"\t"+str(ncol_vals)+"\n")
            z_output.write(header)
            sp_output.write(header)
            continue
        lineList = line.split("\t")

        vals = []
        in_vals =[]
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


    for outline in z_outlines:
        z_output.write(outline + "\n")


    input_file.close()


    #importing data as a matrix
    def strip_first_col(fname, delimiter=None):
        with open(fname, 'r') as fin:
            for line in fin:
                try:
                    yield line.split(delimiter, 1)[1]
                except IndexError:
                    continue

    data = np.loadtxt(strip_first_col(options.input), skiprows=1)

    sp_matrix = (stats.spearmanr(data))


    n=1
    for line in sp_matrix[0]:
        sp_line = np.ndarray.tolist(line)
        header_list = header.split("\t")
        sp_line.insert(0,header_list[n])
        for i in sp_line:
            sp_output.write(str(i).replace("\n", "") + "\t")
        sp_output.write("\n")
        n = n + 1

    z_output.close()
    sp_output.close()
    sys.exit(0)


#############
# FUNCTIONS #
#############

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line


#################
# END FUNCTIONS #
#################
if __name__ == "__main__":
    main()
