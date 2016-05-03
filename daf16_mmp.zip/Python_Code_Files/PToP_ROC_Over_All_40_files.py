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

# Takes in two lists rankedList and geneAndAllRankingsList as well as 
# the name of the KO gene and the "testCount" for that gene, which is
# incremented for each d value used. 
def rankForAllD(rankedList, geneAndAllRankingsList, koGene, testCount):
	num = len(rankedList)

	i = 0
	while (i < num):
		# If the gene is the KO gene and is not in the geneAndAllRankingsList
		# add the gene name to the list and then add the genes rank for this
		# value of d to the list. The rank for the gene will be i + 1 as the
		# rankedList is in acsending order starting from rank 1.
		if ((rankedList[i][0] == koGene) and (rankedList[i][0] not in geneAndAllRankingsList)):
			testCount = testCount + 1
			geneAndAllRankingsList.append(rankedList[i][0])
			geneAndAllRankingsList.insert(testCount,i+1) # insert at element testCount, the value i + 1
		# If the gene is the KO gene and is in the geneAndAllRankingsList
		# simply add the genes rank for the value of d to the list.
		elif (rankedList[i][0] == koGene):
			testCount = testCount + 1
			geneAndAllRankingsList.insert(testCount,i+1) # insert at element testCount, the value i + 1
		i = i + 1
	
	# This list, will be updated each time a new value of d is used
	# Starting with the rank for the gene at d = 0. The list will 
	# have the form [GeneName, rank for d = 0, rank for d = 0.05... 
	# rank for d = 1] testCount is the number of values of d that
	# have been used. The geneAndAllRankingsList only contains the ranking
	# for the KO gene at each value of d.
	return geneAndAllRankingsList, testCount

# rankedList has the form of [geneID, ranking value, gene name, geneID...]
# and koGene is the knocked out gene the experiment was based on. This method gives a
# roc score based on the ranking value of each gene compared against a predicted value.
def createValidityScores(rankedList, koGene, validityTrueList, validityScoreList):
	num = len(rankedList)
	i = 0
	while (i < num):
		# for the knocked out gene, the predicted value should be 1, as
		# it was knocked out, all other genes are given a predicted value
		# of 0
		if (rankedList[i][0] == koGene):
			validityTrueList.append(1)
		else:
			validityTrueList.append(0)
			
		tempValue = (rankedList[i][1])
		validityScoreList.append(tempValue)
		
		i = i + 1
	
	return validityTrueList, validityScoreList

# Write the header information to a specified file, then close the file
# for safety. This appends to file and does not override the files contents.	
def writeHeaderToFile(outputFile):
	with open(outputFile, "a") as this_file:
		this_string = ["GeneName d=0 d=0.05 d=0.1 d=0.15 d=0.2 d=0.25 d=0.3 d=0.35 d=0.4 d=0.45 d=0.5 d=0.55 d=0.6 d=0.65 d=0.7 d=0.75 d=0.8 d=0.85 d=0.9 d=0.95 d=1 \n"]
		this_file.writelines(this_string)
	this_file.close()	

# Takes in a list, 	geneAndAllRankingsList in the form of
# [geneName, Rank1, Rank2, Rank3... Rank20] where rank1 to rank20 are 
# the genes ranking for each value of d measured. As well as the desired
# output file, a count of d values measured (testCount), and a check value 
# (countPerGene)
# Writes the ranking values for all values of d, for each gene		
def writeResultsToFile(geneAndAllRankingsList, outputFile, testCount, CountPerGene, fileCount):
	num2 = len(geneAndAllRankingsList)
	this_ranking = ""
	this_string = ""
	
	with open(outputFile, "a") as this_file:
		# Write out the ranking data for each gene
		i = 1
		while (i < num2):
			if ((testCount == 21) and (CountPerGene == 0)):
				this_ranking = [str(geneAndAllRankingsList), "\n"]
				this_file.writelines(this_ranking)
				CountPerGene = CountPerGene + 1
			i = i + 1
	
		this_file.close()
		
		return CountPerGene, fileCount	
		
# writes the roc score to the file.		
def writeRocScore(validityTrueList, validityScoreList, outputFile):
	y_true = np.array(validityTrueList) 
	y_score = np.array(validityScoreList)
	rocScore = roc_auc_score(y_true, y_score)
	
	with open(outputFile, "a") as this_file:		
		this_string = ["\n", str(rocScore), " "]
		this_file.writelines(this_string)	
	this_file.close()

# Takes in the desired file to read data from, the name of the ko gene,
# the desired output file and the parameter d. The methods needed to 
# run the algorithm and output the data are then run in order as required.	
def runRanking(connectionFilePath, expressionFilePath, koGene, graphFile, outputFile, d, geneAndAllRankingsList, testCount, CountPerGene, fileCount, validityTrueList, validityScoreList):
	fullGeneList, gene1List, gene2List, connectionList = readConnectionFile(connectionFilePath);
	G = makeGraph(fullGeneList);
	G = connectGraph(G, gene1List, gene2List, connectionList, graphFile);
	exprDataList = readExpressionFile(expressionFilePath, fullGeneList);
	normExprDataList = getNormalisedExprData(exprDataList);
	rankingValueList = geneRank(fullGeneList, exprDataList, normExprDataList, G, d);
	rankedList = sortByRanking(fullGeneList, rankingValueList);
	geneAndAllRankingsList, testCount = rankForAllD(rankedList, geneAndAllRankingsList, koGene, testCount);
	#CountPerGene, fileCount = writeResultsToFile(geneAndAllRankingsList, outputFile, testCount, CountPerGene, fileCount);
	validityTrueList, validityScoreList = createValidityScores(rankedList, koGene, validityTrueList, validityScoreList);
	
	return geneAndAllRankingsList, testCount, CountPerGene, fileCount, validityTrueList, validityScoreList

# Takes in the desired file to read data from, the index for GO ID's in
# the given file, the KO gene name, the desired graph output file name,
# the desired ranking output file name, a list containing current rankings
# for the gene, a value showing the number of d values measured, and a check value
# Runs runRanking for all values of d for the given gene.		
def main(connectionFilePath, expressionFilePath, koGene, graphFile, outputFile, geneAndAllRankingsList, testCount, CountPerGene, fileCount):
	#connectionFilePath = input('Give the full file path for the connection file: ')
	#expressionFilePath = input('Give the full file path for the expression file: ')
	#koGene = input('Give the name of the KO gene: ')
	#graphFile = input('What do you want to call the graph file? include .xml ')
	#outputFile = input('What do you want to call the output file for the ranking data? ')
	#runRanking(connectionFilePath, expressionFilePath, koGene, graphFile, outputFile, 0.95);
	
	i = 0.0
	while (i <= 1.02):
		geneAndAllRankingsList, testCount, CountPerGene, fileCount = runRanking(connectionFilePath, expressionFilePath, koGene, graphFile, outputFile, i, geneAndAllRankingsList, testCount, CountPerGene, fileCount);
		i = i + 0.05
	
	print (geneAndAllRankingsList)
	return fileCount
	
#writeHeaderToFile("PToPRankings.txt");	

# Hard-coded values for graph file and output file names as well
# as a list of all KO genes to be tested. Also hard-coded is the
# index for the GO ID's in the given input files. This is done to
# quickly run the ranking for all 40 genes for all values of d for 
# testing and validation purposes and is not designed for end users.
fileCount = 0
testCount = 0
CountPerGene = 0
geneAndAllRankingsList = []
validityTrueList = []
validityScoreList = []
geneConnectionPath = "C:\ThirdYear\Dissertation\daf16_mmp.zip\PtoPDataFiles\\"
geneExprPath = "C:\ThirdYear\Dissertation\daf16_mmp.zip\AmmendedDataFiles\\"
testGeneList = ["Abca1", "Btk", "Cav1", "Cav3", "Cftr", "Clcn1", "Cnr1", "Emd", "Epas1", "Esrra", "Gap43", "Gnmt", "Hdac1", "Hdac2", "Hsf4", "Hspa1a", "Il6", "Lhx1", "Lhx8", "Lmna", "Mbnl1", "Mst1r", "Myd88", "Nos3", "Phgdh", "Pmp22", "Ppara", "Pthlh", "Rab3a", "Rasgrf1", "Rbm15", "Runx2", "Scd1", "Slc26a4", "Srf", "Tcf7", "Tgm2", "Zc3h12a", "Zfp36", "Zfx"]
endOfConnectionPath = "_100nn_stringNet.csv"
endOfExprPath = "_100nn.tsv"
graphPath = "_PToP_graph.xml"	
	
i = 0
num = 40
j = 0.0
while (j < 1.02):
	validityTrueList = []
	validityScoreList = []
	while (i < num):
		geneAndAllRankingsList, testCount, CountPerGene, fileCount, validityTrueList, validityScoreList = runRanking(geneConnectionPath + testGeneList[i] + endOfConnectionPath, geneExprPath + testGeneList[i] + endOfExprPath, testGeneList[i], testGeneList[i] + graphPath, "Rankings.txt", j, geneAndAllRankingsList, testCount, CountPerGene, fileCount, validityTrueList, validityScoreList);
		fileCount = 0
		geneAndAllRankingsList = []
		testCount = 0
		CountPerGene = 0
		i = i + 1

	writeRocScore(validityTrueList, validityScoreList, "Rankings.txt");	
	i = 0
	j = j + 0.05