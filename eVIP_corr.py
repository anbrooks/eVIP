
import sys
import optparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
from scipy import stats

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
            print ("%s option not supplied") % option
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
    opt_parser.add_option("--exp_table",
                          dest="exp_table",
                          type="string",
                          help="""file with filtered and normalized gene expression data table""",
                          default=None)

    opt_parser.add_option("--output",
                          dest="output",
                          type="string",
                          help="""name of output file """,
                          default=None)


    (options, args) = opt_parser.parse_args()




    #validate command line arguments
    opt_parser.check_required("exp_table")
    opt_parser.check_required("output")

    exp_table_file = open(options.exp_table)
    output_file = open(options.output + ".csv","w")

############
# END_MAIN #
############

matrix = np.loadtxt(open("test.csv","rb"),delimiter=",",skiprows=1, usecols=[1,2,3,4])
print(matrix)

print(stats.zscore(matrix))
