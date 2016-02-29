import networkx as nx

# https://networkx.github.io/documentation/latest/reference/classes.graph.html

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
				#print(goList)   proof that this part works
				goList = []	# reset the GO terms list to null 
			else:  
				goList.append(go)

connectedList = []
while (j < len(geneList)):	
	if (i != j):
		#print("This part works") # proof that the idea works
		#print(geneList[i])
		#print(geneList[j])
		if(bool(set(geneList[i]) & set(geneList[j]))):
			if(geneList[i-1] not in connectedList):
				#print("hello")
				connectedList.append(geneList[i-1])
				G.add_edge(geneList[i-1],geneList[j-1])
			elif(geneList[j-1] not in connectedList): 
				#print("goodbye")
				connectedList.append(geneList[j-1])
				G.add_edge(geneList[i-1],geneList[j-1])
			else:
				print("yolo")			
				G.add_edge(geneList[i-1],geneList[j-1])
				#print(geneList[i-1] + ", " + geneList[j-1])
			
			#print("success")  proof that this part works
			i = i + 2
			j = j + 2
	i = i + 2
	j = j + 2
	
nx.write_graphml(G, "testGraph.xml")	
#print(G.edges())      #proof that the gene connections work
	
#print(geneList)	 Shows the list creates properly with all genes	
#print("done") # to show me it is finished loading NOW REDUNDANT