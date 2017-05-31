import mygene
import urllib2
import argparse
import json

#this script gets a current list of pathways available in kegg
#gets the official gene names of each gene in each pathway
#and also converts the official gene names into ensembl ids

##########
#  MAIN  #
##########
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-out", required=True, help="name of output json file that will contain kegg pathway data taken")
    args = parser.parse_args()

    unique_list = []
    all_kegg_pathways = []
    gene_dict_list = []
    pway_gene_dict ={}

    #opening kegg sie with pathway info
    response = urllib2.urlopen('http://rest.kegg.jp/link/hsa/pathway').read()

    for item in response.strip().split("\n"):
        line = item.split()
        if line[0] in unique_list:
            pass
        else:
            unique_list.append(line[0])

    #reformatting the pathway names
    for item in unique_list:
        new_item= item.replace("path:", "")
        all_kegg_pathways.append(new_item)


    #creating dict with each pathway and its genes in ensembl id

    # short_pways = ['hsa01100', 'hsa04141', 'hsa00010']

    for pway in all_kegg_pathways:
        pway_gene_dict[pway] = kegg_to_ensembl_list(pway)

    print pway_gene_dict

    #writing to file
    with open(args.out, "wb") as out_file:
        json.dump(pway_gene_dict, out_file)



############
# END MAIN #
############


#############
# FUNCTIONS #
#############

def kegg_to_ensembl_list(hsa_num):
    #this function takes a kegg pathway number, finds the genes present in the pathway,
    #and converts the genes to ensembl ids

    hsa_num = str(hsa_num.strip('"'))

    #saving kegg website data
    response = urllib2.urlopen('http://rest.kegg.jp/get/'+(hsa_num.strip('"'))).read()

    ##parsing the kegg site

    #creating tuple for start and end indices of the genes
    gene_idx_range = tuple()


    #finding the indices

    check_points = ('COMPOUND','REFERENCE','KO_PATHWAY')

    for index, item in enumerate(response.split()):
        if item == 'GENE':
            gene_idx_range += (index,)

        if item in check_points:
            gene_idx_range += (index,)
            break

    print "\n"
    print "Pathway:"
    print hsa_num

    if len(gene_idx_range) == 1:
        print "No genes in pathway."
        return


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

#################
# END FUNCTIONS #
#################

if __name__ == "__main__": main()
