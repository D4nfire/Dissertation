import networkx as nx

#Parsing data from the gene annotation file to extract a list of
#genes with their GO 

#Reading and parsing the file
geneList = []
tempList = []
goList = []
count = 0
i = 1
j = 3

G = nx.Graph()
with open('C:\ThirdYear\Dissertation\Other helpful documents\Gene_ontology_annotations.txt') as infile:
	for line in infile:
		if (line.startswith('SGD')):
			tempList = line.split("	")
			gene = tempList[2]
			go = tempList[4]					#THIS WORKS!
			#tempGene = (gene + " " + go)
			#print(tempGene) # proof it works
			
			# if not in genelist add new gene to list, then
			# add the go to the gene (can be a list itself
			# else, gene already exists, match and add go
			# to that gene
			# use these genes (within the genelist) as nodes
			
			
			# If the gene isn't in the list AND it is the first gene
			# then do the folowing
			if ((gene not in geneList) and (count == 0)): 
				count = count + 1  # so that it runs for the first gene only
				goList.append(go)	# add go term to a list
				geneList.append(gene) # add the gene to a different list
				G.add_node(gene) 	# make the gene a node in graph
			elif (gene not in geneList):  # runs for subsequent genes
				geneList.append(goList)		# add the GO terms list to the genelist
				geneList.append(gene)	# add the current gene to the list
				G.add_node(gene)	# make the current gene another node in the graph
				goList = []	# reset the GO terms list to null
			#else:  maybe something else here

while (j < len(geneList)):	
	if (i != j):
		#print("This part works") # proof that the idea works
		if (bool(set(geneList[i]) & set(geneList[j]))):
			print("success") # this never seems to be called so far
			i = i + 2
			j = j + 2
	i = i + 2
	j = j + 2
	
#print(geneList)	 Shows the list creates properly with all genes	
print("done") # to show me it is finished loading