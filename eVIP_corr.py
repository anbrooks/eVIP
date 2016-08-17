
import sys
import optparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
from scipy import stats
import cmap.io.gct as gct
import cmap.io.plategrp as grp
import gp


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
    #    if getattr(self.values, option.dest) is None:
    #        print ("%s option not supplied") % option
    #        self.print_help()
    #        sys.exit(1)

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

    # opt_parser.add_option("--output",
    #                       dest="output",
    #                       type="string",
    #                       help="""name of output file """,
    #                       default=None)


    (options, args) = opt_parser.parse_args()

    #validate command line arguments
    #opt_parser.check_required("exp_table")
    #opt_parser.check_required("output")

    #exp_table_file= open(options.exp_table)
    #output_file = open(options.output + ".gct","w")

    #turning expression table into matrix
    #matrix = np.loadtxt(open("test.csv","rb"),delimiter=",",skiprows=1, usecols=[])




    matrix = np.loadtxt(open("test.csv","rb"),delimiter=",",skiprows=1, usecols=[1,2,3,4])




    expmatrix_to_zscore(matrix)

    print(matrix.shape)
    print(matrix)

    #exp_table.close()
    sys.exit(0)




############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def expmatrix_to_zscore(matrix):
    print(matrix)
    #getting zscore matrix
    zscore_mat = stats.zscore(matrix)
    print(zscore_mat)

    #spearman rank pairwise correlatin
    spearm_mat = stats.spearmanr(zscore_mat)
    print(spearm_mat)

def create_gct(output):
    


#################
# END FUNCTIONS #
#################
if __name__ == "__main__": main()
