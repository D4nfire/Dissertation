import random
import numpy as np
import networkx as nx
import scipy as sp
from sklearn.metrics import roc_auc_score
from sklearn import metrics
#import ExpressionData 		#ATM this will run the script as is, before doing this one

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.loadmat.html # to read matlab file
# https://networkx.github.io/documentation/latest/reference/classes.graph.html
# https://networkx.github.io/documentation/latest/reference/generated/networkx.linalg.graphmatrix.adjacency_matrix.html
# http://pythoncentral.io/how-to-sort-a-list-tuple-or-object-with-sorted-in-python/
# http://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html

#Parsing data from the gene annotation file to extract a list of
#genes with their GO 
def readFile():
	#Reading and parsing the file
	geneList = []
	tempList = []
	goList = []
	count = 0

	with open('C:\ThirdYear\Dissertation\Other_helpful_documents\Gene_ontology_annotations.txt') as infile:
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
				# If the gene isn't in the list then add the previous go List
				# to the geneList. Reset the go list and add first go element
				# for this new gene. Add gene to both lists as above.
				elif (gene not in geneList):  # runs for subsequent genes
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
	return geneList

# add each gene to another list for later use
def createList(geneList):
	gene_ex_r = []
	i = 0
	while (i < len(geneList)):
		gene_ex_r.append(geneList[i])
		i = i + 2
	return gene_ex_r
	
def makeGraph(geneList):
	G = nx.Graph()
	i = 0
	while (i < len(geneList)):
		G.add_node(geneList[i])		# add a node for each gene
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
				#print("This part works") # proof that the idea works
				#print(geneList[i])
				#print(geneList[j])
				if(bool(set(geneList[i]) & set(geneList[j]))):			
					G.add_edge(geneList[i-1],geneList[j-1])
					#print(geneList[i-1] + ", " + geneList[j-1])		#Proof it works
			
					#print("success")  proof that this line is reached
				j = j + 2
			j = j + 2
		i = i + 2 	# THE FINAL GENE HAS NO CONNECTIONS, FIX IT
		k = k + 2 # file dropped size a fair bit.
		j = k

	#nx.write_graphml(G, "anotherTestGraph4.xml")	 # Writes the graph to a file
	
	#This will print all the genes and the number of links they have in the graph.
	#i = 0
	#while (i < len(geneList)):
	#	print("The degree for gene" , geneList[i] , "is:" , G.degree(geneList[i]))
	#	i = i + 2
	
	# which can be imported to cytoscape to visualise the graph.
	A = nx.adjacency_matrix(G) 	# Creates an adjacency matrix of the graph above.
	#print (A)		#proof it works
	#print(G.edges())      #proof that the gene connections work	
	#print(geneList)	 #Shows the list creates properly with all genes
	return G
	#print("done") # to show me this part is finished	
				
# This is done purely to get some expression value in
# to do the rest of the algorithm, the ex values here
# are randomly generated floats from -4 to +4. This
# is put into a list in the form of gene, ex, gene, ex etc
def AddEx(gene_ex_r):
	i = 0
	numGenes = len(gene_ex_r) # The number of genes
	num = numGenes*2
	while (i < num):
		temp_ex = random.uniform(-4.0, 4.0)
		temp_ex = abs(temp_ex)
		gene_ex_r.insert(i+1 , temp_ex)
		temp_ex = 0 # A precaution
		i = i + 2
	return gene_ex_r
	#print("done") # to show me this part is finished

#	This finds the sum of all expression values
def sumOfEx(gene_ex_r):
	num = len(gene_ex_r)
	temp_ex = 0
	sumOfEx = 0
	i = 0
	while (i < num):
		temp_ex = abs(gene_ex_r[i+1])
		sumOfEx = sumOfEx + temp_ex
		temp_ex = 0 # A precaution
		i = i + 2
	return sumOfEx
	
# Create the initial ranking of all genes (temp_r)
# Initialises the ranking of each gene, making the
# list have the form of gene, ex, r, gene, ex, r, etc...
# Then add a sumOfConnection element for later use
# making the list have the form of gene, ex, r, sumOf,
# gene, ex, r, sumOf, etc...
def addInitialRanking(gene_ex_r, sumOfEx):
	numGenes = len(gene_ex_r) # The number of genes
	num2 = numGenes*2
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
	return gene_ex_r
	#print("done") # to show me this part is finished

def addSumOf(gene_ex_r):
	sumOfConnection = 0 # saves time later
	numGenes = len(gene_ex_r)
	num2 = numGenes*1.25
	i = 0
	while (i < num2):
		gene_ex_r.insert(i+4 , sumOfConnection)
		i = i + 5
		
	#print (gene_ex_r) 	# gene, expression value, normalised ex value, rank, sumOf, etc...
	return gene_ex_r
	#print("done") # to show me this part is finished

# Run num3 iterations, updating every genes rank
# per iteration
def geneRank(gene_ex_r):
	rankComparisonList = []
	tempRankList = []
	num3 = len(gene_ex_r)
	num4 = len(gene_ex_r)
	d = 0.8
	connectionValue = 0
	i = 0 
	j = 0
	k = 0
	rank = 0
	count = 0
	converged = False
	while ((i < num3) and (converged != True)):
	#while (converged != True):
		# if all genes are compared to all others and it still
		# hasn't converged then start again with the first gene
		#if (i >= len(gene_ex_r)):
		#	i = 0
		#if (count == 0):
		#	rankComparisonList = []
		#	count = count + 1
		#else:
		#	rankComparisonList.append(tempRankList)
		#	tempRankList = [] # resets the list
		while (j < num4):
			if (i != j):
				# If the two genes are connected, hasEdge = 1
				if (G.has_edge(gene_ex_r[i],gene_ex_r[j])):
					hasEdge = 1
				else:
					hasEdge = 0
				# temp_connection is the wij rj[n-1] / degi part, for this gene 
				# The if statement is done so that should a gene have no connections
				# for any reason, the program won't break. The results however may
				# not be as expected.
				if (G.degree(gene_ex_r[i]) == 0):
					temp_connection = (hasEdge * (gene_ex_r[i+1])) / 1 
				else:
					temp_connection = (hasEdge * (gene_ex_r[i+1])) / (G.degree(gene_ex_r[i])) 
				# sumOfConnection is the current sum of the above part, for i iterations
				sumOfConnection = (gene_ex_r[j+4]) + temp_connection
				gene_ex_r[j+4] = sumOfConnection  # update the above sum value for this gene
				connectionValue = (1-d)*(gene_ex_r[i+2]) # this forms the (1-d)exj part   # this takes the normalised ex value
				rank = connectionValue + (d * gene_ex_r[j+3]) 		# this is the ranking of gene j after								
				gene_ex_r[j+3] = rank	# update the ranking		# i iterations
				connectionValue = 0  # Resets the value as a precaution
				tempRankList.append(rank)
				if (count <= 1): 
					converged = False
				else:
					if(bool(set(rankComparisonList[k]) & set(rankComparisonList[k-1]))):
						converged = True
			j = j + 5
		i = i + 5
		j = 0
		k = k + 1
	#print (gene_ex_r) 	# gene, expression value, normalised ex value, rank, sumOf, etc...
	print ("The algorithm converged:" , converged , "after" , k , "itterations")
	return gene_ex_r

def getThreshold():
	i = 0
	threshold = 0
	while (i < len(gene_ex_r)):
		threshold = threshold + gene_ex_r[i+3]
		i = i + 5
		
	j = i/5
	threshold = threshold / j
	return threshold
	
# Creates a list for each gene, it's associated ranking value, the 
# normalised expression value for the gene and a prediction of 1 or 0
# based on it's ranking score, used for later analysis. This then
# adds this to a list containing all genes with their rankings. This list 
# is in the form of [[Gene,RankingValue,normExValue,predicted],
# [Gene,RankingValue,normExValue,predicted]...]
# each gene is then printed with it's rank and ranking value. The list
# of all genes and ranking values is returned for later use.
def sortAndPrintRanking(gene_ex_r, threshold):
	num5 = len(gene_ex_r)	
	rankedList = []
	rankPerGeneList = [] 
	i = 0
	while (i < num5):
			rankPerGeneList.append(gene_ex_r[i]) # This is the gene name
			rankPerGeneList.append(gene_ex_r[i+3]) # This is the ranking score
			rankPerGeneList.append(gene_ex_r[i+1])	# This is the normalised expression value
			# If the ranking score is above a threshold then add an
			# element to the list containing 1, otherwise add 0, this
			# is per gene and is for later use by ROC.
			if (gene_ex_r[i+3] > threshold): 
				rankPerGeneList.append(1)
			else:
				rankPerGeneList.append(0)
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
	return rankedList

def testValidity(rankedList):
	y_true = []
	y_score = []
	y_prediction = []
	num6 = 6381
	i = 0
	while (i < num6):
		tempTrueValue = (rankedList[i][1]) # This is the ranking score for gene i
		tempPredValue = (rankedList[i][2]) # This is the normalised expression value for gene i
		tempPrediction = (rankedList[i][3]) # This is the prediction for gene i
		y_true.append(tempTrueValue)
		y_score.append(tempPredValue)
		y_prediction.append(tempPrediction)
		i = i + 1

	y_true = np.array(y_true)
	#print (y_true)
	#print (y_prediction)
	y_score = np.array(y_score)
	y_prediction = np.array(y_prediction)
	print ("The raw auc value is:" , metrics.auc(y_true, y_score))
	print("The roc auc value is: {:0.3f}".format(roc_auc_score(y_prediction, y_true)))
		
geneList = readFile();
gene_ex_r = createList(geneList);
G = makeGraph(geneList);
G = connectGraph(G, geneList);
gene_ex_r = AddEx(gene_ex_r);
sumOfEx = sumOfEx(gene_ex_r);
gene_ex_r = addInitialRanking(gene_ex_r, sumOfEx);
gene_ex_r = addSumOf(gene_ex_r);
gene_ex_r = geneRank(gene_ex_r);
threshold = getThreshold();
rankedList = sortAndPrintRanking(gene_ex_r, threshold);
testValidity(rankedList);