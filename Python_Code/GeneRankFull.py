import random
import numpy
import networkx as nx
import scipy as sp
#import ExpressionData 		#ATM this will run the script as is, before doing this one

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.loadmat.html # to read matlab file
# https://networkx.github.io/documentation/latest/reference/classes.graph.html
# https://networkx.github.io/documentation/latest/reference/generated/networkx.linalg.graphmatrix.adjacency_matrix.html
# http://pythoncentral.io/how-to-sort-a-list-tuple-or-object-with-sorted-in-python/

#Parsing data from the gene annotation file to extract a list of
#genes with their GO 

#Reading and parsing the file
geneList = []
tempList = []
goList = []
gene_ex_r = []
count = 0
i = 1
j = 1
k = 1

G = nx.Graph()
with open('C:\ThirdYear\Dissertation\Other helpful documents\Gene_ontology_annotations.txt') as infile:
	for line in infile:
		if (line.startswith('SGD')):
			tempList = line.split("	")
			gene = tempList[2]
			go = tempList[4]					
			#tempGene = (gene + " " + go)
			#print(tempGene) # proof it works
			
			# if not in genelist add new gene to list, then
			# add the go to the gene (can be a list itself
			# else, gene already exists, match and add go
			# to that gene
			# use these genes (within the genelist) as nodes
			
			
			# If the gene isn't in the list AND it is the first gene
			# then do the folowing. This adds the go term to the goList, 
			# and the gene to the GeneList, as well as gene_ex_r list
			# for later use. 
			if ((gene not in geneList) and (count == 0)): 
				count = count + 1  # so that it runs for the first gene only
				goList.append(go)	# add go term to a list
				geneList.append(gene) # add the gene to a different list
				gene_ex_r.append(gene)
				G.add_node(gene) 	# make the gene a node in graph
			# If the gene isn't in the list then add the previous go List
			# to the geneList. Reset the go list and add first go element
			# for this new gene. Add gene to both lists as above.
			elif (gene not in geneList):  # runs for subsequent genes
				geneList.append(goList)		# add the GO terms list to the genelist
				geneList.append(gene)	# add the current gene to the list
				gene_ex_r.append(gene)
				G.add_node(gene)	# make the current gene another node in the graph
				goList = []	# reset the GO terms list to null 
				goList.append(go)
			# If the gene is in the list then simply add the go term
			# to the goList for that gene
			else:  
				goList.append(go)

while (i < len(geneList)):	
	while (j < len(geneList)):
		if (i != j):
			#print("This part works") # proof that the idea works
			#print(geneList[i])
			#print(geneList[j])
			if(bool(set(geneList[i]) & set(geneList[j]))):			
				G.add_edge(geneList[i-1],geneList[j-1])
				#print(geneList[i-1] + ", " + geneList[j-1])		#Proof it works
			
				#print("success")  proof that this line is reached
			j = j + 2
		j = j + 2
	i = i + 2
	k = k + 2 # file dropped size a fair bit.
	j = k
	
nx.write_graphml(G, "anotherTestGraph4.xml")	 # Writes the graph to a file
# which can be imported to cytoscape to visualise the graph.
A = nx.adjacency_matrix(G) 	# Creates an adjacency matrix of the graph above.
#print (A)		#proof it works
#print(G.edges())      #proof that the gene connections work
	
#print(geneList)	 #Shows the list creates properly with all genes	
#print("done") # to show me this part is finished				
				
# This is done purely to get some expression value in
# to do the rest of the algorithm, the ex values here
# are randomly generated floats from -4 to +4. This
# is put into a list in the form of gene, ex, gene, ex etc
i = 0
numGenes = len(gene_ex_r) # The number of genes
num = numGenes*2
temp_ex = 0
sumOfEx = 0
while (i < num):
	temp_ex = random.uniform(-4.0, 4.0)
	temp_ex = abs(temp_ex)
	gene_ex_r.insert(i+1 , temp_ex)
	sumOfEx = sumOfEx + temp_ex
	temp_ex = 0 # A precaution
	i = i + 2
	
#print("done") # to show me this part is finished

# Create the initial ranking of all genes (temp_r)
# Initialises the ranking of each gene, making the
# list have the form of gene, ex, r, gene, ex, r, etc...
# Then add a sumOfConnection element for later use
# making the list have the form of gene, ex, r, sumOf,
# gene, ex, r, sumOf, etc...
sumOfConnection = 0 # saves time later
num2 = numGenes*4
temp_r = 0
i = 0
sumOfInitialGenes = 0
while (i < num2): # Ranking.py has the alternative temp_r method
	temp_r = gene_ex_r[i+1]/sumOfEx # GeneRankLoad.py has the other alternative 
	gene_ex_r.insert(i+2 , temp_r) # This will form the initial ranking, the 1-norm value for later use
	gene_ex_r.insert(i+3 , temp_r) # This will be the actual ranking
	sumOfInitialGenes = sumOfInitialGenes + temp_r
	temp_r = 0		# resets the value to 0, ready for the next gene
	i = i + 4
#print (sumOfInitialGenes) # Gives roughly 1.000000000024 as an answer
#print (gene_ex_r) 	# gene, expression value, normalised ex value, rank, etc...
#print("done") # to show me this part is finished

num2 = numGenes*5
i = 0
while (i < num2):
	gene_ex_r.insert(i+4 , sumOfConnection)
	i = i + 5

#print (gene_ex_r) 	# gene, expression value, normalised ex value, rank, sumOf, etc...
#print("done") # to show me this part is finished

# Run num3 iterations, updating every genes rank
# per iteration
num3 = 2000
num4 = len(gene_ex_r)
d = 0.5
connectionValue = 0
i = 0 
j = 0
rank = 0
while (i < num3):
	while (j < num4):
		if (i != j):
			# If the two genes are connected, hasEdge = 1
			if (G.has_edge(gene_ex_r[i],gene_ex_r[j])):
				hasEdge = 1
			else:
				hasEdge = 0
			# temp_connection is the wij rj[n-1] / degi part, for this gene
			temp_connection = (hasEdge * (gene_ex_r[i+1])) / (G.degree(gene_ex_r[i]))
			# sumOfConnection is the current sum of the above part, for i iterations
			sumOfConnection = (gene_ex_r[j+4]) + temp_connection
			gene_ex_r[j+4] = sumOfConnection  # update the above sum value for this gene
			connectionValue = (1-d)*(gene_ex_r[i+2]) # this forms the (1-d)exj part   # this takes the normalised ex value
			rank = connectionValue + (d * gene_ex_r[j+3]) 		# this is the ranking of gene j after								
			gene_ex_r[j+3] = rank	# update the ranking		# i iterations
			connectionValue = 0  # Resets the value as a precaution
		j = j + 5
	i = i + 5
	j = 0
#print (gene_ex_r) 	# gene, expression value, normalised ex value, rank, sumOf, etc...
	
num5 = numGenes*5	
rankedList = []
rankPerGeneList = []
i = 0
while (i < num5):
		rankPerGeneList.append(gene_ex_r[i])
		rankPerGeneList.append(gene_ex_r[i+3])
		rankedList.append(rankPerGeneList)
		rankPerGeneList = []
		i = i + 5
		
def getKey(item):
	return item[1]
rankedList = sorted(rankedList, key = getKey, reverse = True)

num6 = len(rankedList)
i = 0
while (i < num6):
	print(rankedList[i][0] , " is ranked: " , i+1, " with a ranking value of: ", rankedList[i][1])
	i = i + 1

