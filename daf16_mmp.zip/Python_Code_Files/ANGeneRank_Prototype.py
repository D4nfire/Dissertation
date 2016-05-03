import random
import networkx as nx

# web sources for some parts of the code. 		

# Using Networkx package:
# http://networkx.readthedocs.io/en/networkx-1.11/tutorial/tutorial.html 
# Sorting a list in python used in sortAndPrintRanking:
# http://pythoncentral.io/how-to-sort-a-list-tuple-or-object-with-sorted-in-python/ 
# GeneRank is based on the algorithm developed by Morrison et al. in the 2005 paper
# "GeneRank: Using search engine technology for the analysis of microarray experiments" 

#Parsing data from the gene annotation file to extract a list of
#genes with their GO annotations
def readFile():
	geneList = [] # a list of genes
	tempList = [] # a temporary list for the elements of each line
	goList = [] # a list of GO annotations
	count = 0
	
	# Each gene and the go annotations associated with that gene are added
	# to the geneList in order of reading.
	with open('C:\ThirdYear\Dissertation\daf16_mmp.zip\DataFiles\Gene_ontology_annotations.txt') as infile:
		for line in infile:
			# if the line has gene data
			if (line.startswith('SGD')): 
				tempList = line.split("	")
				gene = tempList[2] # gene name
				go = tempList[4] # GO annotations	
			
				# If the gene isn't in the list AND it is the first gene
				# then add the go term to the goList, and the gene to the GeneList
				if ((gene not in geneList) and (count == 0)): 
					count = count + 1  # increment so that it runs for the first gene only
					goList.append(go)	
					geneList.append(gene) 
				# Otherwise, if the gene isn't in the list then add the previous go List
				# to the geneList, now in the order of [gene name, go terms]. Add the next
				# gene to the list and Reset the go list before adding to it for the next gene.
				elif (gene not in geneList):  
					geneList.append(goList)		# add the GO terms list to the genelist
					geneList.append(gene)	# add the current gene to the list
					goList = []	# reset the GO terms list to null 
					goList.append(go)
				# If the gene is in the list then simply add the go term
				# to the goList for that gene
				else:  
					goList.append(go)
			# This is done so that the last gene still has it's GO terms added
			elif (line == '\n'):	
				geneList.append(goList)
	# geneList takes the form of [gene name, [go terms], gene name, [go terms]...]
	return geneList
	
# Create the network graph with each gene name as a node
def makeGraph(geneList):
	G = nx.Graph()
	i = 0
	while (i < len(geneList)):
		G.add_node(geneList[i])	# add a node for each gene
		i = i + 2
	return G
	
#Takes in graph and geneList and adds the gene connections to the graph
def connectGraph(G, geneList):
	i = 1
	j = 1
	k = 1
	while (i < len(geneList)):	
		while (j < len(geneList)):
			if (i != j):
				# if the set of go terms for gene1 and gene2 have a match
				# then add an edge between these two genes
				if(bool(set(geneList[i]) & set(geneList[j]))):			
					G.add_edge(geneList[i-1],geneList[j-1])
				j = j + 2
			j = j + 2
		i = i + 2 	

	nx.write_graphml(G, "anotherTestGraph4.xml")	 # Writes the graph to a file
	return G
	
# This is done purely to get some expression value in
# for the rest of the algorithm to run, the ex values here
# are randomly generated floats from -4 to +4. This
# is put into a list in the form of gene, ex, gene, ex etc
def RandGenerateEx(geneList):
	gene_ex_r = []
	numGenes = len(geneList) # The number of genes
	num = numGenes

	# add all the gene names to gene_ex_r
	i = 0
	while (i < num):
		gene_ex_r.append(geneList[i]) # add the gene name to the list
		
		i = i + 2

	# add expr values to gene_ex_r
	i = 0
	while (i < num):
		temp_ex = random.uniform(-4.0, 4.0)
		temp_ex = abs(temp_ex) # testing that the abs function works as expected
		gene_ex_r.insert(i+1, temp_ex) # add the expression value to the list
		temp_ex = 0 # A precaution
		
		i = i + 2
		
	# gene_ex_r has the form of [gene name, expr value, gene name... etc]
	return gene_ex_r

# calculate the sum of all expression values.
def sumOfEx(gene_ex_r):
	num = len(gene_ex_r)
	temp_ex = 0
	sumOfExpr = 0
	
	i = 1
	while (i < num):
		temp_ex = abs(gene_ex_r[i])
		sumOfExpr = sumOfExpr + temp_ex
		temp_ex = 0 # A precaution
		i = i + 2	
	return sumOfExpr
	
# add the initial ranking and an element for the actual ranking for each gene. 
def addInitialRanking(gene_ex_r, sumOfEx):
	numGenes = len(gene_ex_r) # The number of genes
	num2 = numGenes*2 # list doubles in size again
	temp_r = 0
	
	i = 1
	while (i < num2): 
		temp_r = gene_ex_r[i]/sumOfEx # calculates the normalised expression value
		gene_ex_r.insert(i+1 , temp_r) # This will form the initial ranking, the 1-norm value of the gene
		gene_ex_r.insert(i+2 , temp_r) # This will form the actual ranking score to be updated
		temp_r = 0		# resets the value to 0, ready for the next gene
		i = i + 4
	
	# gene_ex_r list takes the form of [gene name, expression value, normalised expression value, 
	# ranking score, gene name etc...]
	return gene_ex_r

# add an element for use within the GeneRank algorithm
def addSumOfConnection(gene_ex_r):
	sumOfConnection = 0 
	numGenes = len(gene_ex_r)
	num2 = numGenes*1.25 # list increases in size by 25%
	
	i = 0
	while (i < num2):
		gene_ex_r.insert(i+4 , sumOfConnection)
		i = i + 5
		 	
	# gene_ex_r has the form [gene name, expr, normalised expr, ranking score, sumOf, gene name, etc...]
	return gene_ex_r

# The generank algorithm iterates over every genes, for each gene
# the ranking of all other genes are updated.
def geneRank(gene_ex_r, G):
	num3 = len(gene_ex_r)
	d = 0.65
	connectionValue = 0
	rank = 0
	
	i = 0 
	j = 0
	while (i < num3):
		while (j < num3):
			if (i != j):
				# If the two genes are connected, hasEdge = 1
				if (G.has_edge(gene_ex_r[i],gene_ex_r[j])):
					hasEdge = 1
				else:
					hasEdge = 0
				# temp_connection is the wij rj[n-1] / degi part, for this gene 
				# The if statement is done so that should a gene have no connections
				# for any reason, the program won't break.
				if (G.degree(gene_ex_r[i]) == 0):
					temp_connection = (hasEdge * (gene_ex_r[i+1])) / 1 
				else:
					temp_connection = (hasEdge * (gene_ex_r[i+1])) / (G.degree(gene_ex_r[i])) 
				# sumOfConnection holds the the current sum of the above part, for i iterations, for each gene
				sumOfConnection = (gene_ex_r[j+4]) + temp_connection
				gene_ex_r[j+4] = sumOfConnection  # update the above sum value for this gene
				connectionValue = (1-d)*(gene_ex_r[i+2]) # this forms the (1-d)exj part   # this takes the normalised ex value
				rank = connectionValue + (d * gene_ex_r[j+3]) 		# this is the ranking of gene j after i iterations							
				gene_ex_r[j+3] = rank	# update the ranking score for this gene	
				
				connectionValue = 0  # Resets the value as a precaution
				rank = 0
				
			j = j + 5
		i = i + 5
		j = 0
	
	# gene_ex_r has the form of [gene name, expr, normalised expr, ranking score, sumOf, gene name, etc...]
	return gene_ex_r

# creates a sorted list of genes and their rankings.
def sortAndPrintRanking(gene_ex_r):
	num4 = len(gene_ex_r)	
	rankedList = []
	rankPerGeneList = [] # a temporary list used for each gene in turn
	
	i = 0
	while (i < num4):
			rankPerGeneList.append(gene_ex_r[i]) # gene name
			rankPerGeneList.append(gene_ex_r[i+3]) # genes ranking score
			rankedList.append(rankPerGeneList)
			rankPerGeneList = []
			i = i + 5
	
	# sorts the list based on ranking score
	def getKey(item):
		return item[1]
	rankedList = sorted(rankedList, key = getKey, reverse = True)

	# prints the genes in rank order with the associated ranking score
	num5 = len(rankedList)
	i = 0
	while (i < num5):
		print(rankedList[i][0] , " is ranked: " , i+1, " with a ranking value of: ", rankedList[i][1])
		i = i + 1

# main method which calls all other methods
def main():	
	geneList = readFile();
	G = makeGraph(geneList);
	G = connectGraph(G, geneList);
	gene_ex_r = RandGenerateEx(geneList);
	sumOfExpr = sumOfEx(gene_ex_r);
	gene_ex_r = addInitialRanking(gene_ex_r, sumOfExpr);
	gene_ex_r = addSumOfConnection(gene_ex_r);
	gene_ex_r = geneRank(gene_ex_r, G);
	sortAndPrintRanking(gene_ex_r);
	
main();