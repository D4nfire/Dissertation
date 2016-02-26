#Parsing data from the gene annotation file to extract a list of
#genes with their GO 

#Reading and parsing the file
geneAnnotationList=[]
tempList=[]
with open('C:\ThirdYear\Dissertation\Other helpful documents\Gene_ontology_annotations.txt') as infile:
	for line in infile:
		v1, v2, v3, v4, v5 = line.split('\t')
		print(v3 + " " + v5)