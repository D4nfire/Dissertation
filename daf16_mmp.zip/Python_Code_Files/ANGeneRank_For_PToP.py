import networkx as nx
from sklearn.metrics import roc_auc_score

# web sources for some parts of the code. 		

# Using Networkx package:
# http://networkx.readthedocs.io/en/networkx-1.11/tutorial/tutorial.html 
# Sorting a list in python used in sortAndPrintRanking:
# http://pythoncentral.io/how-to-sort-a-list-tuple-or-object-with-sorted-in-python/
# Creating the roc score
# http://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html
# GeneRank is based on the algorithm developed by Morrison et al. in the 2005 paper
# "GeneRank: Using search engine technology for the analysis of microarray experiments" 

# Takes in the protein-protein connection file
def readConnectionFile(filepath):
	#Reading and parsing the file
	fullGeneList = []
	gene1List = []
	gene2List = []
	connectionList = []

	with open(filepath) as infile:
		next(infile)
		for line in infile:
			tempList = line.split(",")
			# if the line doesn't contain the gene names then skip that line.
			if ((tempList[0] == "") or (tempList[1] == "")):
				next(infile)
			# otherwise parse the data to fullGeneList
			elif ((tempList[0] not in fullGeneList) and (tempList[1] not in fullGeneList)):
				fullGeneList.append(tempList[0]) # Gene name
				fullGeneList.append(tempList[1]) # Gene name
				
			elif (tempList[0] not in fullGeneList):
				fullGeneList.append(tempList[0]) # Gene name
				
			elif (tempList[1] not in fullGeneList):
				fullGeneList.append(tempList[1]) # Gene name

			gene1List.append(tempList[0]) # for later use
			gene2List.append(tempList[1])
			
			# If the two genes are connected by the desired source
			# "experiments" then add a value of 1 to the list. 
			# Otherwise add 0. The list will be in the order of reading
			if (tempList[7] == 0):
				connectionList.append(0)
			else:
				connectionList.append(1)
	
	return fullGeneList, gene1List, gene2List, connectionList

# For each gene name add a node to the graph G, return the graph for later us
# fullGeneDataList contains all the gene names 
def makeGraph(fullGeneList):
	G = nx.Graph()
	
	i = 0
	while (i < len(fullGeneList)):
		G.add_node(fullGeneList[i])	# add a node for each gene name
		i = i + 4
		
	return G

# For each gene pair, if there is a connection then add 
# an edge in the graph between these genes	
def connectGraph(G, gene1List, gene2List, connectionList, graphFile):
	num = len(connectionList)
	
	i = 0
	while (i < num):
		# if the two genes corresponding connection value is 1
		if (connectionList[i] == 1):
			G.add_edge(gene1List[i],gene2List[i])
		i = i + 1

	#Write the graph to a file, which can be read by Cytoscape and shown visually.
	#nx.write_graphml(G, graphFile) 
	
	return G	

# Read in and parse data from a given file. This file must be in a certain format 
# where expr data must be the 5th element
def readExpressionFile(filePath, fullGeneList):
	#Reading and parsing the file
	exprDataList = []
	tempList = []
	num = len(fullGeneList)

	# initialise the expression values to 0
	i = 0
	while (i < num):
		exprDataList.append(0)
		i = i + 1
	
	with open(filePath) as infile:
		next(infile)
		for line in infile:
			tempList = line.split("	")
			# if the line doesn't contain the gene name or expr value then skip that line.
			if ((tempList[6] == "") or (tempList[4] == "") or (tempList[6] not in fullGeneList)):
				next(infile)
			# otherwise parse the data to exprDataList
			else:
				# put the expr data in the same gene order as fullgeneList
				i = 0
				while (i < num):
					if (fullGeneList[i] == tempList[6]):
						exprDataList[i] = (abs(float(tempList[4])))
					i = i + 1
	
	# this list contains the expression data for each gene, in the same
	# order as the fullgeneList
	return exprDataList
	
# calculates the normalised expr value for each gene
def getNormalisedExprData(exprDataList):
	normExprDataList = []
	sumOfEx = 0 
	temp_norm_ex = 0
	num = len(exprDataList)
	
	i = 0
	while (i < num):
		sumOfEx = sumOfEx + exprDataList[i]
		i =  i + 1 
	
	j = 0
	while (j < num):
		temp_norm_ex = exprDataList[j] / sumOfEx
		normExprDataList.append(temp_norm_ex)
		j = j + 1
		
	return normExprDataList

# The main algorithm takes in three lists, geneIDList, exprDataList and normExprDataList
# as well as the Graph G and parameter d. The algorithm is run such that each gene is
# given an initial ranking score, the normalised expression value of that gene. The ranking
# score is updated each iteration. The algorithm runs for n iterations where n is the number 
# of genes. !The data in the three lists must be in the same gene order, such that 
# geneIDList[i], exprDataList[i] and normExprDataList[i] all belong to the same gene!
def geneRank(fullGeneList, exprDataList, normExprDataList, G, d):
	sumOfConnectionList = []
	rankingValueList = []
	num = len(fullGeneList)

	i = 0 
	j = 0
	# The algorithm runs for "num" iterations, with num being the number of genes.
	while (i < num):
		# each gene has it's ranking updated each iteration
		while (j < num):
			# On the first iteration add some data for each gene
			if (i == 0):
				sumOfConnectionList.append(0) # A value for each gene which is updated after each iteration
				rankingValueList.append(normExprDataList[j]) # initialise ranking with the normalised expr value
			# If two genes are the same then they can't be compared and updated so skip to the next gene
			elif (i != j): 
				# If the two genes are connected in graph G, hasEdge = 1
				if (G.has_edge(fullGeneList[i],fullGeneList[j])): 
					hasEdge = 1
				else:
					hasEdge = 0
				# temp_connection is a value based on (wij * r(j)^[i-1]) / outDegree(i) or 
				# 1/0 based on connection, multiplied by the rank of gene j in the previous iteration,
				# all divided by the outDegree of gene i.
				if (G.degree(fullGeneList[i]) == 0):
					temp_connection = (hasEdge * (rankingValueList[i-1])) / 1 
				else:
					temp_connection = (hasEdge * (rankingValueList[i-1])) / (G.degree(fullGeneList[i])) 
				# sumOfConnection[j] is the current sum of the above part for gene j, after i iterations
				sumOfConnectionList[j] = sumOfConnectionList[j] + temp_connection
				connectionValue = (1-d)*exprDataList[j]
				rankValue = connectionValue + d*(sumOfConnectionList[j]) # The final ranking value
				rankingValueList[j] = rankValue # update the ranking score for gene j
				
			j = j + 1
		i = i + 1
		j = 0 # reset j so that every gene is compared again next iteration

	return rankingValueList # List of all genes ranking values in the order they were put in

# Takes in two lists, fullGeneList and rankingValueList and then
# combines the data into one list, rankedList which has the form of
# [gene name, ranking value, gene name...] This new list is then
# sorted based on the ranking value, so the first gene in the list will
# be the one with the highest ranking value.
def sortByRanking(fullGeneList, rankingValueList):
	geneAndRankList = []
	rankedList = []
	num = len(fullGeneList)
	
	i = 0
	while (i < num):
		geneAndRankList.append(fullGeneList[i])
		geneAndRankList.append(rankingValueList[i])
		rankedList.append(geneAndRankList) 
		geneAndRankList = [] # reset list to null for the next genes data
		i = i + 1
	
	# sorts the list based on ranking value.
	def getKey(item):
		return item[1]
	rankedList = sorted(rankedList, key = getKey, reverse = True)
	
	return rankedList

# rankingValueList has the form of [gene name, ranking value, gene name...]
# koGene is the knocked out gene the experiment was based on. This method gives a
# roc score based on the ranking value of each gene compared against a predicted value.
def testValidity(rankedList, koGene):
	y_true = [] # A list for the true values given by the algorithm
	y_score = [] # A list for the predicted values
	num6 = len(rankedList)	
	
	i = 0
	while (i < num6):
		# for the knocked out gene, the predicted value should be 1, as
		# it was knocked out, all other genes are given a predicted value
		# of 0
		if (rankedList[i][0] == koGene):
			y_true.append(1)
		else:
			y_true.append(0)
			
		tempValue = (rankedList[i][1])
		y_score.append(tempValue)
			
		i = i + 1
		
	y_true = np.array(y_true) 
	y_score = np.array(y_score)
	rocScore = roc_auc_score(y_true, y_score)
	
	return rocScore	
	
# Takes in the final ranking list, the measured roc Score, the parameter d
# and the desired output file name. Then Writes desired output to a file, this
# is a header, followed by ranking information for each gene, followed by the
# roc score.
def writeResultsToFile(rankedList, rocScore, d, outputFile):
	num2 = len(rankedList)
	this_ranking = ""
	this_string = ""
	
	with open(outputFile, "a") as this_file:
	
		this_string = ["Ranking and roc information for variable d = ", str(d), "\n\n"]
		this_file.writelines(this_string)
	
		# Write out the ranking data for each gene
		i = 0
		while (i < num2):
			this_ranking = [str(rankedList[i][0]), " is ranked: ", str(i+1) , " with a ranking value of: ", str(rankedList[i][1]), "\n"]
			this_file.writelines(this_ranking)
			i = i + 1
	
		this_string = ["\nThe roc score is: " , str(rocScore), "\n\n\n\n"]
		this_file.writelines(this_string)
		this_file.close()	

# Takes in the desired file to read data from, the name of the ko gene,
# the desired output file and the parameter d. The methods needed to 
# run the algorithm and output the data are then run in order as required.
def runRanking(connectionFilePath, expressionFilePath, koGene, graphFile, outputFile, d):
	fullGeneList, gene1List, gene2List, connectionList = readConnectionFile(connectionFilePath);
	G = makeGraph(fullGeneList);
	G = connectGraph(G, gene1List, gene2List, connectionList, graphFile);
	exprDataList = readExpressionFile(expressionFilePath, fullGeneList);
	normExprDataList = getNormalisedExprData(exprDataList);
	rankingValueList = geneRank(fullGeneList, exprDataList, normExprDataList, G, d);
	rankedList = sortByRanking(fullGeneList, rankingValueList);
	rocScore = testValidity(rankedList, koGene);
	writeResultsToFile(rankedList, rocScore, d, outputFile);

# Runs the file and asks for the file path information and the KO gene
# name. d is set to 0.65 but can be cahnged.	
def main():
	connectionFilePath = input('Give the full file path for the connection file: ')
	expressionFilePath = input('Give the full file path for the expression file: ')
	koGene = input('Give the name of the KO gene: ')
	graphFile = input('What do you want to call the graph file? include .xml ')
	outputFile = input('What do you want to call the output file for the ranking data? ')
	runRanking(connectionFilePath, expressionFilePath, koGene, graphFile, outputFile, 0.65);
	
main();