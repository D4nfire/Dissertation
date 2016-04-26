import random
import numpy as np
import networkx as nx
import scipy as sp
import csv
from sklearn.metrics import roc_auc_score

# References for code	

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.loadmat.html # to read matlab file
# https://networkx.github.io/documentation/latest/reference/classes.graph.html
# https://networkx.github.io/documentation/latest/reference/generated/networkx.linalg.graphmatrix.adjacency_matrix.html
# http://pythoncentral.io/how-to-sort-a-list-tuple-or-object-with-sorted-in-python/
# http://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html

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
				fullGeneList.append(tempList[0]) # Gene ID
				fullGeneList.append(tempList[1]) # Gene ID
				
			elif (tempList[0] not in fullGeneList):
				fullGeneList.append(tempList[0]) # Gene ID
				
			elif (tempList[1] not in fullGeneList):
				fullGeneList.append(tempList[1]) # Gene ID

			gene1List.append(tempList[0])
			gene2List.append(tempList[1])
			
			if (tempList[7] == 0):
				connectionList.append(0)
			else:
				connectionList.append(1)
	
	return fullGeneList, gene1List, gene2List, connectionList

# For each gene ID add a node to the graph G, return the graph for later us
# fullGeneDataList has the form of [Gene ID, GO ID list, expr data, gene name, Gene ID...] 
def makeGraph(fullGeneList):
	G = nx.Graph()
	
	i = 0
	while (i < len(fullGeneList)):
		G.add_node(fullGeneList[i])	# add a node for each gene
		i = i + 4
		
	return G

# For each gene, check the genes list of GO ID's against every other gene
# if two genes both share one or more GO ID's then add a link between these
# genes in the graph G.		
	
def connectGraph(G, gene1List, gene2List, connectionList, graphFile):
	num = len(connectionList)
	
	i = 0
	while (i < num):
		if (connectionList[i] == 1):
			G.add_edge(gene1List[i],gene2List[i])
		i = i + 1

	#Write the graph to a file, which can be read by Cytoscape and shown visually.
	#nx.write_graphml(G, graphFile) 
	
	return G	

# Read in and parse data from a given file. This file must be in a certain format 
# where geneID must be the first element of each line, expr data must be the 5th element,
# gene name must be the 7th element and the GO id's must be the 14th, 15th and 16th
# elements. This information is put into a list called fullGeneDataList.
#genes with their GO 
def readExpressionFile(filePath, fullGeneList):
	#Reading and parsing the file
	exprDataList = []
	tempList = []
	num = len(fullGeneList)

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
			# otherwise parse the data to fullGeneDataList
			else:
				# put the expr data in the same gene order as fullgeneList
				i = 0
				while (i < num):
					if (fullGeneList[i] == tempList[6]):
						exprDataList[i] = (abs(float(tempList[4])))
					i = i + 1
	
	# This list has the form of [Gene ID, GO ID list, expr data, gene name, Gene ID...] 
	return exprDataList
	
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

# The main algorithm takes in three lists, fullGeneList, exprDataList and normExprDataList
# as well as the Graph G and parameter d. Using this information the algorithm is run such
# that each gene is given a ranking based on a ranking score or value. The ranking value,
# and therefore the overall ranking of each gene is updated each iteration. The algorithm
# runs for n iterations where n is the number of genes. !The data in the three lists must be
# in the same gene order, such that fullGeneList[i], exprDataList[i] and normExprDataList[i]
# all belong to the same gene!
def geneRank(fullGeneList, exprDataList, normExprDataList, G, d):
	sumOfConnectionList = []
	rankingValueList = []
	num3 = len(fullGeneList)

	i = 0 
	j = 0
	# The algorithm runs num3 iterations, with num3 being the number of genes.
	while (i < num3):
		# each gene has it's ranking updated each iteration
		while (j < num3):
			# On the first iteration add some data for each gene
			if (i == 0):
				sumOfConnectionList.append(0) # A value for each gene which is updated after each iteration
				rankingValueList.append(normExprDataList[j]) # initialise ranking with the normalised expr value
			# If two genes are the same then they can't be compared and updated so skip
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
				###temp_connection = (hasEdge * (rankingValueList[i-1])) / (G.degree(geneIDList[i]))
				# sumOfConnection is the current sum of the above part for gene j, after i iterations
				sumOfConnectionList[j] = sumOfConnectionList[j] + temp_connection
				connectionValue = (1-d)*normExprDataList[j]
				rankValue = connectionValue + d*(sumOfConnectionList[j]) # The final ranking value
				rankingValueList[j] = rankValue # List of all genes ranking values
				
			j = j + 1
		i = i + 1
		j = 0 # reset j so that every gene is compared again next iteration

	return rankingValueList # List of all genes ranking values	

# Takes in three lists, geneIDList, rankingValueList and geneNameList, 
# combines the data into one list, rankedList which has the form of
# [geneID, ranking value, gene name, geneID...] This new list is then
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

def rankForAllD(rankedList, geneAndAllRankingsList, koGene, testCount):
	num = len(rankedList)
	i = 0
	while (i < num):
		if ((rankedList[i][0] == koGene) and (rankedList[i][0] not in geneAndAllRankingsList)):
			testCount = testCount + 1
			geneAndAllRankingsList.append(rankedList[i][0])
			geneAndAllRankingsList.insert(testCount,i+1)
		elif (rankedList[i][0] == koGene):
			testCount = testCount + 1
			geneAndAllRankingsList.insert(testCount,i+1)
		i = i + 1
	return geneAndAllRankingsList, testCount

def writeHeaderToFile(outputFile):
	with open(outputFile, "a") as this_file:
		this_string = ["GeneName d=0 d=0.05 d=0.1 d=0.15 d=0.2 d=0.25 d=0.3 d=0.35 d=0.4 d=0.45 d=0.5 d=0.55 d=0.6 d=0.65 d=0.7 d=0.75 d=0.8 d=0.85 d=0.9 d=0.95 d=1 \n"]
		this_file.writelines(this_string)
	this_file.close()	
	
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
	
def runRanking(connectionFilePath, expressionFilePath, koGene, graphFile, outputFile, d, geneAndAllRankingsList, testCount, CountPerGene, fileCount):
	fullGeneList, gene1List, gene2List, connectionList = readConnectionFile(connectionFilePath);
	G = makeGraph(fullGeneList);
	G = connectGraph(G, gene1List, gene2List, connectionList, graphFile);
	exprDataList = readExpressionFile(expressionFilePath, fullGeneList);
	normExprDataList = getNormalisedExprData(exprDataList);
	rankingValueList = geneRank(fullGeneList, exprDataList, normExprDataList, G, d);
	rankedList = sortByRanking(fullGeneList, rankingValueList);
	geneAndAllRankingsList, testCount = rankForAllD(rankedList, geneAndAllRankingsList, koGene, testCount);
	CountPerGene, fileCount = writeResultsToFile(geneAndAllRankingsList, outputFile, testCount, CountPerGene, fileCount);
	
	return geneAndAllRankingsList, testCount, CountPerGene, fileCount
	
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
	
writeHeaderToFile("PToPRankings.txt");	
	
fileCount = 0
geneAndAllRankingsList = []
validityTrueList = []
validityScoreList = []
geneConnectionPath = "C:\ThirdYear\Dissertation\Python_Code\GR_V1.9\PtoPDataFiles\\"
geneExprPath = "C:\ThirdYear\Dissertation\Python_Code\GR_V1.9\DataFiles\\"
testGeneList = ["Abca1", "Btk", "Cav1", "Cav3", "Cftr", "Clcn1", "Cnr1", "Emd", "Epas1", "Esrra", "Gap43", "Gnmt", "Hdac1", "Hdac2", "Hsf4", "Hspa1a", "Il6", "Lhx1", "Lhx8", "Lmna", "Mbnl1", "Mst1r", "Myd88", "Nos3", "Phgdh", "Pmp22", "Ppara", "Pthlh", "Rab3a", "Rasgrf1", "Rbm15", "Runx2", "Scd1", "Slc26a4", "Srf", "Tcf7", "Tgm2", "Zc3h12a", "Zfp36", "Zfx"]
endOfConnectionPath = "_100nn_stringNet.csv"
endOfExprPath = "_100nn.tsv"
graphPath = "_PToP_graph.xml"	
	
num = 40
i = 0
while (i < num):
	geneAndAllRankingsList = []
	testCount = 0
	CountPerGene = 0
	fileCount = main(geneConnectionPath + testGeneList[i] + endOfConnectionPath, geneExprPath + testGeneList[i] + endOfExprPath, testGeneList[i], testGeneList[i] + graphPath, "PToPRankings.txt", geneAndAllRankingsList, testCount, CountPerGene, fileCount);
	i = i + 1