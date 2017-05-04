import argparse

########
# MAIN #
########
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gage_same", help="GAGE results when ran with same.dir=TRUE; represents the up or down regulated pathways")
    parser.add_argument("-gage_diff", help="GAGE results when ran with same.dir=FALSE; represents the dysregulated pathways")
    parser.add_argument("-eVIPP_summary", help="eVIPP_summary.txt file from run_eVIP with -eVIPP option")
    parser.add_argument("-out", help="Name of output file which will contain the eVIPP calls and the GAGE results")
    args = parser.parse_args()

#getting list of hsa pathway numbers from eVIPP summary output
    header = None
    hsa_list = []
    dict = {}
    for line in open(args.eVIPP_summary, "U"):
        if line.startswith("pathway"):
            header = line
            continue
        else:
            hsa_list.append(line.split()[0])
            dict[line.split()[0]]= line.split()[1:]

    new_hsa_list=[]
#removing duplicate hsa numbers from the list
    for i in hsa_list:
        if i not in new_hsa_list:
            new_hsa_list.append(i)

#finding the values in the gage_same file
    for hsa in new_hsa_list:
        for line in open(args.gage_same, "U"):
            val = line.split(",")[0]
            if val.split()[0][1:] == hsa:
                dict[hsa] += ' '.join(val.split()[1:])[:-1],line.split(",")[4],line.split(",")[10]

        for line in open(args.gage_diff, "U"):
            val = line.split(",")[0]
            if val.split()[0][1:] == hsa:
                dict[hsa] += (line.split(",")[4],)

    with open(args.out, "w+") as out:
        out.writelines(header.strip("\n") + ("\t")+ "pathway_name" +("\t") + "gage_same_greater_qval" + ("\t")+ "gage_same_less_qval" + ("\t")+ "gage_diff_dys_qval" +("\n"))
        for key in dict:
            out.write(key + "\t" + "\t".join(dict[key]) + "\n")





############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatLine(line):
    line = line.replace("\r", "")
    return line

#################
# END FUNCTIONS #
#################

if __name__ == "__main__": main()
