import csv

with open("C:\ThirdYear\Dissertation\Other helpful documents\Gene_ontology_annotations.txt") as tsv:
	for line in csv.reader(tsv, dialect="excel-tab"):
		Gene = line[2]
		GO = line[4]
		tempGene = (Gene + " " + GO)
		print(tempGene)