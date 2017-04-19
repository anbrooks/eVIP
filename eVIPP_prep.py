import argparse
import urllib2
import sys
import mygene
import os

########
# MAIN #
########
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-hsa_num", help="kegg pathway number, for running eVIP on a single pathway")
    parser.add_argument("-hsa_list", help="file containing a list of kegg pathway numbers, for running eVIP on several pathways at once ")
    parser.add_argument("-data_file", help="file with data to extract pathway genes from")
    parser.add_argument("-o", help = "name of output directory to create")

    args = parser.parse_args()

    hsa_num = args.hsa_num
    hsa_list = open(args.hsa_list, "r")
    data = open(args.data_file, "r")

    #make eVIP output directory
    out_dir = args.o
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if hsa_num:
        #run on the one pathway
        kegg_to_data(hsa_num)

    if hsa_list:
        pathway_list = []
        for line in hsa_list:
            pathway_list.append(line.strip())

        for pathway in pathway_list:
            out_file = open(out_dir + "/" +pathway +".txt", "w")
            ensembl_list = kegg_to_ensembl_list(pathway)
            gene_extraction(ensembl_list, open(args.data_file), out_file)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############

def kegg_to_ensembl_list(hsa_num):
    #this function takes a kegg pathway number, finds the genes present in the pathway,
    #and converts the genes to ensembl ids

    #saving kegg website data
    response = urllib2.urlopen('http://rest.kegg.jp/get/'+hsa_num).read()

    ##parsing the kegg site

    #creating tuple for start and end indices of the genes
    gene_idx_range = tuple()


    #finding the indices

    check_points = ('COMPOUND','REFERENCE')

    for index, item in enumerate(response.split()):
        if item == 'GENE':
            gene_idx_range += (index,)

        if item in check_points:
            gene_idx_range += (index,)
            break

    print "\n"
    print "Pathway:"
    print hsa_num

    #segmenting the data based off the start and end indices
    kegg_items = (response.split()[gene_idx_range[0]:gene_idx_range[1]])

    #getting the gene names
    gene_names = []
    for item in kegg_items:
        if ";" in item:
            item = item.replace(";", "")
            gene_names.append(item)

    print "Genes in pathway:"
    print gene_names

    #converting from official gene symbol to ensembl id
    mg = mygene.MyGeneInfo()
    gene_dict_list = mg.querymany(gene_names, scopes='symbol',fields='ensembl.gene',species='human')

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



def gene_extraction(ensembl_list, data, output):
# this function extracts the pathway genes from the RNA-seq data
    matched_ids = []
    for line in data:
        if line.startswith("#"):
            output.write(line)
        line = line.split()
        for id in ensembl_list:
            #if id ==line[0]
            if line[0].startswith(id):
                matched_ids.append(id)
                output.write("\t".join(line) + "\n")


    matched_length = str(len(matched_ids))
    ensembl_length = str(len(ensembl_list))
    print matched_length + " of " + ensembl_length + " IDs were found in the data. "



#################
# END FUNCTIONS #
#################

if __name__ == "__main__":
    main()