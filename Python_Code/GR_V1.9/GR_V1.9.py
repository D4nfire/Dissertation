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

# Read in and parse data from a given file. This file must be in a certain format 
# where geneID must be the first element of each line, expr data must be the 5th element,
# gene name must be the 7th element and the GO id's must be the 14th, 15th and 16th
# elements. This information is put into a list called fullGeneDataList.
#genes with their GO 
def readFile(filePath):
	#Reading and parsing the file
	fullGeneDataList = []
	goList = []
	tempList = []

	with open(filePath) as infile:
		next(infile)
		for line in infile:
			tempList = line.split("	")
			# if the line doesn't contain the gene name or expr value then skip that line.
			if ((tempList[6] == "") or (tempList[4] == "")):
				next(infile)
			# otherwise parse the data to fullGeneDataList
			else:
				fullGeneDataList.append(tempList[0]) # Gene ID
				tempGoList1 = (tempList[13].split("///")) #GO ID's
				tempGoList2 = (tempList[14].split("///")) #GO ID's
				tempGoList3 = (tempList[15].split("///")) #GO ID's
				goList = tempGoList1 + tempGoList2 + tempGoList3 # merge the three GO ID lists into one list
				fullGeneDataList.append(goList)
				fullGeneDataList.append(abs(float(tempList[4])))
				fullGeneDataList.append(tempList[6])
				goList = [] # resets the temporary list for the next line
	
	# This list has the form of [Gene ID, GO ID list, expr data, gene name, Gene ID...] 
	#print (fullGeneDataList)
	return fullGeneDataList

# For each gene ID add a node to the graph G, return the graph for later use
# fullGeneDataList has the form of [Gene ID, GO ID list, expr data, gene name, Gene ID...] 
def makeGraph(fullGeneDataList):
	G = nx.Graph()
	
	i = 0
	while (i < len(fullGeneDataList)):
		G.add_node(fullGeneDataList[i])	# add a node for each gene
		i = i + 4
		
	return G

# For each gene, check the genes list of GO ID's against every other gene
# if two genes both share one or more GO ID's then add a link between these
# genes in the graph G.	
def connectGraph(G, fullGeneDataList, graphFile):
	i = 1
	j = 1
	while (i < len(fullGeneDataList)):	
		while (j < len(fullGeneDataList)):
			if (i != j):
				if(bool(set(fullGeneDataList[i]) & set(fullGeneDataList[j]))):			
					G.add_edge(fullGeneDataList[i-1],fullGeneDataList[j-1])
				j = j + 4
			j = j + 4
		i = i + 4
		j = 1

	#Write the graph to a file, which can be read by Cytoscape and shown visually.
	nx.write_graphml(G, graphFile) 
	
	return G

# fullGeneDataList has the form of [geneID, GO ID list, ex_data, geneName, Gene ID...].
# Split the fullGeneDataList into three lists, geneIDList, exprDataList, normExprDataList
# and geneNameList as these are required for the main algorithm and later methods.
def getMultipleLists(fullGeneDataList):
	geneIDList = []
	exprDataList = []
	normExprDataList = []
	geneNameList = []
	sumOfEx = 0 # the sum of all expr values 
	
	i = 0
	j = 0
	while (i<len(fullGeneDataList)):
		geneIDList.append(fullGeneDataList[i])
		exprDataList.append(fullGeneDataList[i+2])
		sumOfEx = sumOfEx + (fullGeneDataList[i+2])
		geneNameList.append(fullGeneDataList[i+3])
		i = i + 4
	
	# normExprDataList contains each genes expr data divided by the sum
	# of all the expr data of all the genes.
	while (j < len(exprDataList)):
		temp_norm_ex = exprDataList[j] / sumOfEx
		normExprDataList.append(temp_norm_ex)
		j = j + 1

	return geneIDList,exprDataList,normExprDataList,geneNameList

# The main algorithm takes in three lists, geneIDList, exprDataList and normExprDataList
# as well as the Graph G and parameter d. Using this information the algorithm is run such
# that each gene is given a ranking based on a ranking score or value. The ranking value,
# and therefore the overall ranking of each gene is updated each iteration. The algorithm
# runs for n iterations where n is the number of genes. !The data in the three lists must be
# in the same gene order, such that geneIDList[i], exprDataList[i] and normExprDataList[i]
# all belong to the same gene!
def geneRank(geneIDList, exprDataList, normExprDataList, G, d):
	sumOfConnectionList = []
	rankingValueList = []
	num3 = len(geneIDList)

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
				if (G.has_edge(geneIDList[i],geneIDList[j])): 
					hasEdge = 1
				else:
					hasEdge = 0
				# temp_connection is a value based on (wij * r(j)^[i-1]) / outDegree(i) or 
				# 1/0 based on connection, multiplied by the rank of gene j in the previous iteration,
				# all divided by the outDegree of gene i.
				if (G.degree(geneIDList[i]) == 0):
					temp_connection = (hasEdge * (rankingValueList[i-1])) / 1 
				else:
					temp_connection = (hasEdge * (rankingValueList[i-1])) / (G.degree(geneIDList[i])) 
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

# fullGeneDataList has the form of [geneID, GO ID list, ex_data, geneName, Gene ID...]
# and koGene is the knocked out gene the experiment was based on. This method gives a
# roc score based on the ranking value of each gene compared against a predicted value.
def testValidity(fullGeneDataList, koGene):
	y_true = [] # A list for the true values given by the algorithm
	y_score = [] # A list for the predicted values
	num6 = len(fullGeneDataList)	
	
	i = 0
	while (i < num6):
		# for the knocked out gene, the predicted value should be 1, as
		# it was knocked out, all other genes are given a predicted value
		# of 0
		if (fullGeneDataList[i+3] == koGene):
			y_true.append(1)
		else:
			y_true.append(0)
			
		tempValue = (fullGeneDataList[i+2])
		y_score.append(tempValue)
			
		i = i + 4
		
	y_true = np.array(y_true) 
	y_score = np.array(y_score)
	rocScore = roc_auc_score(y_true, y_score)
	
	return rocScore

# Takes in three lists, geneIDList, rankingValueList and geneNameList, 
# combines the data into one list, rankedList which has the form of
# [geneID, ranking value, gene name, geneID...] This new list is then
# sorted based on the ranking value, so the first gene in the list will
# be the one with the highest ranking value.
def sortByRanking(geneIDList, rankingValueList, geneNameList):
	geneAndRankList = []
	rankedList = []
	num = len(geneIDList)
	
	i = 0
	while (i < num):
		geneAndRankList.append(geneIDList[i])
		geneAndRankList.append(rankingValueList[i])
		geneAndRankList.append(geneNameList[i])
		rankedList.append(geneAndRankList) 
		geneAndRankList = [] # reset list to null for the next genes data
		i = i + 1
	
	# sorts the list based on ranking value.
	def getKey(item):
		return item[1]
	rankedList = sorted(rankedList, key = getKey, reverse = True)
	
	return rankedList

# Takes in the final ranking list, the measured rocScore, the parameter d
# and the desired output file name. Print to file, the information about
# the next set of data, followed by the data, followed by the rocScore for
# that data. The output file is appended each time and not replaced so multiple
# tests data can be written to the same file.
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
			this_ranking = [str(rankedList[i][0]), " (" , str(rankedList[i][2]) , ") ", "is ranked: ", str(i+1) , " with a ranking value of: ", str(rankedList[i][1]), "\n"]
			this_file.writelines(this_ranking)
			i = i + 1
	
		this_string = ["\nThe roc score is: " , str(rocScore), "\n\n\n\n"]
		this_file.writelines(this_string)
		this_file.close()

# Takes in the desired file to read data from, the name of the ko gene,
# the desired output file and the parameter d. The methods needed to 
# run the algorithm and output the data are then run in order as required.
def runRanking(filepath, koGene, graphFile, outputFile, d):
	#filepath = input('Give the full file path: ')
	#koGene = input('Give the name of the ko gene: ')
	#outputFile = input('What do you want to call the output file? ')
	fullGeneDataList = readFile(filepath);
	G = makeGraph(fullGeneDataList);
	G = connectGraph(G, fullGeneDataList, graphFile);
	geneIDList,exprDataList,normExprDataList,geneNameList = getMultipleLists(fullGeneDataList);
	rankingValueList = geneRank(geneIDList, exprDataList, normExprDataList, G, d);
	rocScore = testValidity(fullGeneDataList, koGene);
	rankedList = sortByRanking(geneIDList, rankingValueList, geneNameList);
	writeResultsToFile(rankedList, rocScore, d, outputFile);

def main():
	# Base tests for each koGene file, this can be removed and the algorithm run
	# purely by calling the geneRank algorithm method with the correct parameters.
	# The output of that would be the ranking value of each gene, in the original
	# order of genes as is in the input file.
	filepath = input('Give the full file path: ')
	koGene = input('Give the name of the ko gene: ')
	graphFile = input('What do you want to call the graph file? include .xml ') 
	outputFile = input('What do you want to call the output file? ')	
	
	i = 0.0
	while (i <= 1.0):
		runRanking(filepath, koGene, graphFile, outputFile, i);
		i = i + 0.05

main();