import sys
import optparse
import os
import pdb
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
    z_output= open(options.z_output, "w+")
    sp_output = open(options.sp_output, "w")

    t_vals=[]
    z_outline = []
    z_outlines = []

    header = None
    for line in input_file:
        line = formatLine(line)
        if line.startswith("#"):
            header = line
            header += "\n"
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

        t_vals = []

    input_file.close()
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