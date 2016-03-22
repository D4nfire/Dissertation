import random
import numpy
import networkx as nx

# https://networkx.github.io/documentation/latest/reference/classes.graph.html
# https://networkx.github.io/documentation/latest/reference/generated/networkx.linalg.graphmatrix.adjacency_matrix.html

# List used to store genes, later expression values,
# rankings and a holding sum value will also be added
# for each gene.
gene_ex_r = []

G = nx.Graph()
with open('C:\ThirdYear\Dissertation\Other helpful documents\Gene_ontology_annotations.txt') as infile:
	for line in infile:
		if (line.startswith('SGD')):
			tempList = line.split("	")
			gene = tempList[2]
			
			gene_ex_r.append(gene)
	
# Read in a pre-existing graph.
G = nx.read_graphml("anotherTestGraph4.xml")	 # reads the graph from a file

A = nx.adjacency_matrix(G) 	# Creates an adjacency matrix of the graph above.
#print (A)		#proof it works
#print(G.edges())      #proof that the gene connections work
	
#print(geneList)	 #Shows the list creates properly with all genes	
print("done") # to show me this part is finished				
				
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
	gene_ex_r.insert(i+1 , temp_ex)
	temp_ex = abs(temp_ex)
	sumOfEx = sumOfEx + temp_ex
	i = i + 2
	
print("done") # to show me this part is finished

# Create the initial ranking of all genes (temp_r)
# Initialises the ranking of each gene, making the
# list have the form of gene, ex, r, gene, ex, r, etc...
# Then add a sumOfConnection element for later use
# making the list have the form of gene, ex, r, sumOf,
# gene, ex, r, sumOf, etc...
sumOfConnection = 0 # saves time later
num2 = numGenes*3
i = 0
while (i < num2): # Ranking.py has the alternative temp_r method
	temp_r = gene_ex_r[i+1]/(numpy.linalg.norm([sumOfEx], ord=1)) 
	gene_ex_r.insert(i+2 , temp_r)
	i = i + 3
#print (gene_ex_r) 	# gene, expression value, rank, etc...
print("done") # to show me this part is finished

num2 = numGenes*4
i = 0
while (i < num2):
	gene_ex_r.insert(i+3 , sumOfConnection)
	i = i + 4

#print (gene_ex_r) 	# gene, expression value, rank, sumOf, etc...
print("done") # to show me this part is finished

# Run num3 iterations, updating every genes rank
# per iteration
num3 = 400
num4 = len(gene_ex_r)
d = 0.5
connectionValue = (1 - d)
i = 0 
j = 0
rank = 0
while (i < num3):
	while (j < num4):
		if (i != j):
			if (G.has_edge(gene_ex_r[i],gene_ex_r[j])):
				hasEdge = 1
			else:
				hasEdge = 0
			temp_connection = (hasEdge * (gene_ex_r[i+1])) / (G.degree(gene_ex_r[i]))
			sumOfConnection = (gene_ex_r[j+3]) + temp_connection
			gene_ex_r[j+3] = sumOfConnection
			rank = connectionValue + (d * gene_ex_r[j+3])
			gene_ex_r[j+2] = rank
		j = j + 4
	i = i + 4
	j = 0

#print (gene_ex_r) 	# gene, expression value, rank, sumOf, etc...	

# Print the ranking value for each gene
#num5 = numGenes*4
#i = 0	
#while (i < num5):
	#print(gene_ex_r[i] , " is ranked " , gene_ex_r[i+2])
	#i = i + 4