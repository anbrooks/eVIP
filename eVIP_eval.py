#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# <Script name>
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb
import csv

#############
# CONSTANTS #
#############
#################
# END CONSTANTS #
#################


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
                          dest="input_table",
                          type="string",
                          help="Input table with mutation impact predictions.",
                          default=None)
#   opt_parser.add_option("-n",
#                         dest="num_predictions",
#                         type="int",
#                         help="Number of predictions in input_table",
#                         default=None)
    opt_parser.add_option("-b",
                          dest="benchmark",
                          type="string",
                          help="Table of benchmark calls of GOF, LOF, or Inert",
                          default=None)
    opt_parser.add_option("--gof_accuracy",
                          dest="gof_accuracy",
                          action="store_true",
                          help="""Indicates that the benchmark is PC9. A separate
                                  accuracy will be reported which allows Inert
                                  predictions to be accurate when the benchmark
                                  is LOF or Inert""",
                          default=False)
    opt_parser.add_option("-o",
                          dest="out_table",
                          type="string",
                          help="""New table which adds on the benchmark
                                  prediction for each allele. The column header will be
                                  the file name""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
#    opt_parser.check_required("-n")
    opt_parser.check_required("-b")
    opt_parser.check_required("-o")

    input_table = open(options.input_table)
    out_table = open(options.out_table, "w")

    benchmark_table = open(options.benchmark)
    gof_accuracy = options.gof_accuracy
    
    mutation2known = parseBenchmark(benchmark_table)

#    num_pred = options.num_predictions
    num_pred = 1
    pred_ctr = [0 for n in range(num_pred)]
    
    if gof_accuracy:
        gof_pred_ctr = [0 for n in range(num_pred)]

    # [{}, {}, {}, {}...num_pred]
    # {} = {known_mut_call:Counts of prediction
    pred_matrix = getPredRecord(num_pred)

    total_calls = [0 for n in range(num_pred)]

    tp_calls = [0 for n in range(num_pred)]
    tn_calls = [0 for n in range(num_pred)]
    fp_calls = [0 for n in range(num_pred)]
    fn_calls = [0 for n in range(num_pred)]

    ctr = 0
    for line in input_table:
        line = formatLine(line)
        lineList = line.split("\t")
        if ctr == 0:
            mut_idx = lineList.index("mut")
            out_table.write(line + "\t" + options.benchmark + "\n")
            ctr += 1
            continue
        
        mut = lineList[mut_idx] 
        if mut not in mutation2known:
            out_table.write(line + "\tNI\n")
            continue
   
        # Write out what the benchmark prediction is
        out_table.write(line + "\t" + mutation2known[mut] + "\n")
 
#        pdb.set_trace()
   
        for i in range(-num_pred, 0):

            if lineList[i] != "NI":
                total_calls[i] += 1        

            if lineList[i] == mutation2known[mut]:
                pred_ctr[i] += 1.0
            else:
                # Allow for COF, DOM-NEG to match with GOF and LOF
                if lineList[i] == "COF" or lineList[i] == "DOM-NEG":
                    if mutation2known[mut] == "GOF" or mutation2known[mut] == "LOF" or mutation2known[mut] == "COF":
                        pred_ctr[i] += 1.0

            if gof_accuracy:
                if mutation2known[mut] == "GOF":
                    if lineList[i] == "GOF" or lineList[i] == "COF":
                        gof_pred_ctr[i] += 1.0
                else:
                    if lineList[i] == "LOF" or lineList[i] == "Inert":
                        gof_pred_ctr[i] += 1.0

            # Record prediction
            pred_matrix[i][mutation2known[mut]][lineList[i]] += 1

            # TP, FP, TN, FN
            if ((mutation2known[mut] == "GOF" or mutation2known[mut] == "LOF" or mutation2known[mut] == "COF") and 
                (lineList[i] == "GOF" or lineList[i] == "LOF" or lineList[i] == "COF" or lineList[i] == "DOM-NEG")):
                tp_calls[i] += 1.0
            if mutation2known[mut] == "Inert" and (lineList[i] == "GOF" or 
                                                   lineList[i] == "COF" or 
                                                   lineList[i] == "LOF" or
                                                   lineList[i] == "DOM-NEG"):
                fp_calls[i] += 1.0
            if mutation2known[mut] == "Inert" and lineList[i] == "Inert":
                tn_calls[i] += 1.0
            if ((mutation2known[mut] == "GOF" or mutation2known[mut] == "LOF" or mutation2known[mut] == "COF") and 
                (lineList[i] == "Inert" or lineList[i] == "NI")):
                fn_calls[i] += 1.0
            
  
    # Print results 
    for i in range(num_pred): 
        try:
            print "Prediction %d accuracy: %.4f" % (i + 1,
                                                    pred_ctr[i]/total_calls[i])
        except:
            print "Prediction %d accuracy: NA" % (i + 1)

    if gof_accuracy:
        print 
        print "GOF accuracy"
        for i in range(num_pred):
            try:
                print "Prediction %d accuracy: %.4f" % (i+1,
                                                    gof_pred_ctr[i]/total_calls[i])
            except:
                print "Prediction %d accuracy: NA" % (i+1)

    # Print sensitivity
    print
    for i in range(num_pred):
        try:
            print "Prediction %d sensitivity: %.4f" % (i+1,
                                                       tp_calls[i]/(tp_calls[i] + fn_calls[i]))
        except:
            print "Prediction %d sensitivity: NA" % (i+1)
            
    print
    # Print specificity
    for i in range(num_pred):
        try:
            print "Prediction %d specificity: %.4f" % (i+1,
                                                       tn_calls[i]/(fp_calls[i] + tn_calls[i]))
        except:
            print "Prediction %d specificity: NA"  % (i+1)

    # Print prediction matrix
    for i in range(num_pred):
        print 
        print "Prediction %d:" % (i + 1)
        print "Known\tGOF_pred\tLOF_pred\tCOF_pred\tDOM-NEG_pred\tInert_pred\tNI_pred"
        print "GOF\t%d\t%d\t%d\t%d\t%d\t%d" % (pred_matrix[i]["GOF"]["GOF"],
                                       pred_matrix[i]["GOF"]["LOF"],
                                       pred_matrix[i]["GOF"]["COF"],
                                       pred_matrix[i]["GOF"]["DOM-NEG"],
                                       pred_matrix[i]["GOF"]["Inert"],
                                       pred_matrix[i]["GOF"]["NI"])
        print "LOF\t%d\t%d\t%d\t%d\t%d\t%d" % (pred_matrix[i]["LOF"]["GOF"],
                                       pred_matrix[i]["LOF"]["LOF"],
                                       pred_matrix[i]["LOF"]["COF"],
                                       pred_matrix[i]["LOF"]["DOM-NEG"],
                                       pred_matrix[i]["LOF"]["Inert"],
                                       pred_matrix[i]["LOF"]["NI"])
        print "COF\t%d\t%d\t%d\t%d\t%d\t%d" % (pred_matrix[i]["COF"]["GOF"],
                                       pred_matrix[i]["COF"]["LOF"],
                                       pred_matrix[i]["COF"]["COF"],
                                       pred_matrix[i]["COF"]["DOM-NEG"],
                                       pred_matrix[i]["COF"]["Inert"],
                                       pred_matrix[i]["COF"]["NI"])
        print "Inert\t%d\t%d\t%d\t%d\t%d\t%d" % (pred_matrix[i]["Inert"]["GOF"],
                                       pred_matrix[i]["Inert"]["LOF"],
                                       pred_matrix[i]["Inert"]["COF"],
                                       pred_matrix[i]["Inert"]["DOM-NEG"],
                                       pred_matrix[i]["Inert"]["Inert"],
                                       pred_matrix[i]["Inert"]["NI"])

    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getPredRecord(num_pred):
    pred_matrix = [{} for n in range(num_pred)]

    for this_dict in pred_matrix:
        this_dict["GOF"] = {"GOF":0,
                            "LOF":0,
                            "COF":0,
                            "DOM-NEG":0,
                            "Inert":0,
                            "NI":0}
        this_dict["LOF"] = {"GOF":0,
                            "LOF":0,
                            "COF":0,
                            "DOM-NEG":0,
                            "Inert":0,
                            "NI":0}
        this_dict["COF"] = {"GOF":0,
                            "LOF":0,
                            "COF":0,
                            "DOM-NEG":0,
                            "Inert":0,
                            "NI":0}
        this_dict["Inert"] = {"GOF":0,
                              "LOF":0,
                              "COF":0,
                              "DOM-NEG":0,
                              "Inert":0,
                              "NI":0}

    return pred_matrix

def parseBenchmark(benchmark_table):
    ctr = 0
    allele2known = {}
    for line in benchmark_table:
        if ctr == 0:
            ctr += 1
            continue
        line = formatLine(line)
        lineList = line.split("\t")
    
        allele2known[lineList[0]] = lineList[1]         
        ctr += 1

    return allele2known

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
