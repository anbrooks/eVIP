import sys
import optparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
from scipy import stats
import pandas as pd

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
    opt_parser.add_option("-i",
                          dest="input",
                          type="string",
                          help="""file with filtered and normalized gene expression data table""",
                          default=None)
    opt_parser.add_option("--z_out",
                          dest="zscore_output",
                          type="string",
                          help="""name of zscore output file """,
                          default=None)
    opt_parser.add_option("--sp_out",
                          dest="spearman_output",
                          type="string",
                          help="""name of spearman output file """,
                          default=None)

    (options, args) = opt_parser.parse_args()

    #validate command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("--z_out")
    opt_parser.check_required("--sp_out")


    input_file = open(options.input)
    z_output_file = open(options.zscore_output + ".csv","w")
    sp_output_file = open(options.spearman_output + ".csv","w")

    #counts # of columns in input file, used for making matrix
    ncols = len(input_file.readline().split(","))

    #TODO define header
    #TODO define genenames

    matrix = np.loadtxt(input_file, delimiter=",",skiprows=0, usecols=range(1,ncols))

    # matrix_file = open('/Users/Alexis/Desktop/matrix.csv', 'w')
    # writer = csv.writer(matrix_file)
    #TODO
    #writer.writerow(header)
    #TODO add genenames

    # for values in matrix:
    #     writer.writerow(values)

    # z_calc = stats.zscore(matrix)
    # z_writer = csv.writer(z_output_file)
    # #z_writer.writerow(header)
    # for z_values in z_calc:
    #     z_writer.writerow(z_values)

    z_writer = csv.writer(z_output_file)
    sp_writer = csv.writer(sp_output_file)

    z_lines = []
    sp_lines = []
    for line in matrix:
        z_line = stats.zscore(line)
        z_lines.append(z_line)
        z_writer.writerow(z_line)
    for z_line in z_lines:
        #sp_line = stats.spearmanr(z_line)
        sp_line = pd.corr(method='spearman')
        sp_lines.append(sp_line)
        sp_writer.writerow(sp_line)



    input_file.close()
    z_output_file.close()
    sp_output_file.close()
#############
# FUNCTIONS #
#############

#################
# END FUNCTIONS #
#################

############
# END_MAIN #
############
if __name__ == "__main__":
    main()
