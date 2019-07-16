import mygene
import argparse
import json

"""
this script takes in a a gene set file from MSigDB
and converts to JSON format to be used with eVIPP
gmt files from:
http://software.broadinstitute.org/gsea/downloads.jsp
"""

##########
#  MAIN  #
##########
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-JSON", required=True, help="location of JSON file to create")
    parser.add_argument("-gmt", help="location of MSigDB gmt file with gene symbols")
    parser.add_argument("-ensembl", action="store_true", help= "Will convert gmt with gene symbols into ensembl ids")
    args = parser.parse_args()

    pway_gene_dict = {}

    with open(args.gmt, "r") as gmt:
        for line in gmt:
            sline = line.split()
            pway = sline[0]
            print pway
            # sline1 is a html link
            ids = sline[2:]

            if args.ensembl:
                #convert ids from gene name to ensembl
                pway_gene_dict[pway] = geneToEnsembl(ids)

            else:
                pway_gene_dict[pway] = ids


    # writing JSON to file
    with open(args.JSON, "wb") as out_file:
        json.dump(pway_gene_dict, out_file)


############
# END MAIN #
############


#############
# FUNCTIONS #
#############

def geneToEnsembl(ids):
    #converting from official gene symbol to ensembl id
    mg = mygene.MyGeneInfo()
    gene_dict_list = mg.querymany(ids, scopes='symbol',fields='ensembl.gene',species='human')

    gene_dicts = []

    for dict in gene_dict_list:
        for key, value in dict.iteritems():
            if key == "ensembl":
                gene_dicts.append(value)

    print "# of gene names:"
    print len(gene_dicts)
    count = 0
    ensembl_list = []
    for d in gene_dicts:
        try:
            value = d.get('gene')
            ensembl_list.append(value)
            count += 1
        except AttributeError:
            for item in d:
                value = item.get('gene')
                ensembl_list.append(value)
                count += 1

    print "# of ensembl ids found:"
    print(count)

    return ensembl_list


#################
# END FUNCTIONS #
#################

if __name__ == "__main__": main()
