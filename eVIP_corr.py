import sys
import optparse
import os
import pdb
import numpy as np
from scipy import stats
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


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
    # opt_parser.add_option("--sp_out",
    #                       dest="sp_output",
    #                       type="string",
    #                       help="""name of output file containing spearman correlation""",
    #                       default=None)

    (options, args) = opt_parser.parse_args()
    opt_parser.check_required("-i")
    opt_parser.check_required("--z_out")
    # opt_parser.check_required("--sp_out")

    input_file = open(options.input)
    # spearman_output = open(options.sp_output, "w")
    z_output_file = open(options.z_output, "w+")


    outlines = []
    mads = []
    #zscores = []
    out_vals =[]
    in_vals =[]
    vals=[]
    t_vals=[]

    header = None
    for line in input_file:
        line = formatLine(line)
        if line.startswith("#"):
            header = line
            z_output_file.write(header + "\n")
            continue

        lineList = line.split()

        for item in lineList[1:]:
            val = float(item)
            in_vals.append(val)
            vals.append(val)
            #print(vals)

        median = robjects.r['mean'](robjects.FloatVector(vals))[0]
        mad = robjects.r['sd'](robjects.FloatVector(vals))[0]



        for v in in_vals:
            print("median")
            print(median)
            print("sd")
            print(mad)
            if mad == 0:
                new_val = (v-median)/CLOSE_TO_ZERO
                t_vals.append(new_val)
                out_vals.append(new_val)
            else:
                new_val = (v-median)/mad
                t_vals.append(new_val)
                out_vals.append(new_val)
            print("newval")
            print(new_val)

        out_vals_str = []
        for v in out_vals:
            out_vals_str.append("%.4f" % v)

        outline = "\t".join(lineList[:1]) + "\t"
        outline += "\t".join(out_vals_str)
        outline += "\n"


        outlines.append(outline)
        #print(outlines)
        mads.append(mad)

    for outline in outlines:
        z_output_file.write(outline)
        #zscores.append(outline)

    #print(lineList)
    #print(lineList[1:])



    input_file.close()
    z_output_file.close()
    #spearman_output.close()
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
