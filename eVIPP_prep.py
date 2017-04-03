import urllib2
import sys
import mygene

hsa_num = sys.argv[1]
data_file = sys.argv[2]
output_file = sys.argv[3]

#saving kegg website data
response = urllib2.urlopen('http://rest.kegg.jp/get/'+hsa_num).read()

#parsing the kegg site

#creating tuple for start and end indices of the genes
gene_idx_range = tuple()

#finding the indices
for index, item in enumerate(response.split()):
    if item == 'GENE' or item == 'COMPOUND':
        gene_idx_range += (index,)

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

#opening file of data
data = open(data_file, "r")

#opening output file
output = open(output_file, "w")

matched_ids = []
for line in data:
    if line.startswith("#"):
        output.write(line)
    line = line.split()
    for id in ensembl_list:
        if id ==line[0]:
            matched_ids.append(id)
            output.write("\t".join(line) + "\n")


matched_length = str(len(matched_ids))
ensembl_length = str(len(ensembl_list))
print matched_length + " of " + ensembl_length + " IDs were found in the data. "
