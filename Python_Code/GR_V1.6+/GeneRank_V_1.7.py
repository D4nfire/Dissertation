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

# Have code started for parsing and using data from new files
# need to sort out the reading of the file, and change some
# parts depending on whether multiple items can be returned or not
# if it does work then add some comments

#Parsing data from the gene annotation file to extract a list of
#genes with their GO 
def readFile():
	#Reading and parsing the file
	fullGeneDataList = []
	goList = []
	tempList = []

	with open('C:\ThirdYear\Dissertation\Python_Code\GR_V1.6\GSE5496.txt') as infile:
		next(infile)
		for line in infile:
			tempList = line.split("	")		
			if ((tempList[6] == "") or (tempList[4] == "")):
				next(infile)
			else:
				fullGeneDataList.append(tempList[0])
				tempGoList1 = (tempList[13].split("///"))
				tempGoList2 = (tempList[14].split("///"))
				tempGoList3 = (tempList[15].split("///"))
				goList = tempGoList1 + tempGoList2 + tempGoList3
				fullGeneDataList.append(goList)
				fullGeneDataList.append(float(tempList[5]))
				fullGeneDataList.append(tempList[6])
				goList = []
	
	#print (fullGeneDataList)
	return fullGeneDataList

def makeGraph(fullGeneDataList):
	G = nx.Graph()
	i = 0
	while (i < len(fullGeneDataList)):
		G.add_node(fullGeneDataList[i])		# add a node for each gene
		i = i + 4
	return G
	
def connectGraph(G, fullGeneDataList):
	i = 1
	j = 1
	k = 1
	while (i < len(fullGeneDataList)):	
		while (j < len(fullGeneDataList)):
			if (i != j):
				if(bool(set(fullGeneDataList[i]) & set(fullGeneDataList[j]))):			
					G.add_edge(fullGeneDataList[i-1],fullGeneDataList[j-1])
				j = j + 4
			j = j + 4
		i = i + 4 	# THE FINAL GENE HAS NO CONNECTIONS, FIX IT
		k = k + 4 # file dropped size a fair bit.
		j = k

	nx.write_graphml(G, "testGraph.xml")
	A = nx.adjacency_matrix(G) 	# Creates an adjacency matrix of the graph above.
	
	return G

# assuming incoming list is in the form of [geneID, go_values, ex_data, geneName].
def getMultipleLists(fullGeneDataList):
	geneIDList = []
	exprDataList = []
	normExprDataList = []
	geneNameList = []
	sumOfEx = 0
	i = 0
	j = 0
	
	while (i<len(fullGeneDataList)):
		geneIDList.append(fullGeneDataList[i])
		exprDataList.append(fullGeneDataList[i+2])
		sumOfEx = sumOfEx + (fullGeneDataList[i+2])
		geneNameList.append(fullGeneDataList[i+3])
		i = i + 4
	
	while (j < len(exprDataList)):
		temp_norm_ex = exprDataList[j] / sumOfEx
		normExprDataList.append(temp_norm_ex)
		j = j + 1
		
	#print (geneIDList)
	#print (exprDataList)
	#print (normExprDataList)
	return geneIDList,exprDataList,normExprDataList,geneNameList

	# CAN USE ADJACENCY MATRIX BUT COULD ALSO USE GRAPH AS geneIDList[i] and geneIDList[j]
	# would be the nodes to check, this is doable. can now call this algorithm with other
	# data more easily.
def geneRank(geneIDList, exprDataList, normExprDataList, G, d):
	rankComparisonList = []
	tempRankList = []
	sumOfConnectionList = []
	rankingValueList = []
	converged = False
	num3 = len(geneIDList)
	num4 = len(geneIDList)
	i = 0 
	j = 0
	k = 0
	count = 0
	
	while (j < num4):
		sumOfConnectionList.append(0) # A value for each gene which is updated after each itteration
		rankingValueList.append(normExprDataList[j]) # initialise ranking with the normalised expr value
		j = j + 1
	j = 0
	
	#while ((i < num3) and (converged != True)):
	while ((i < 100) and (converged != True)):
		# if all genes are compared to all others and it still
		# hasn't converged then start again with the first gene
		if (i >= len(geneIDList)):
			numOfItterations = numOfItterations + 1
			i = 0
		else:
			if (count == 0):
				rankComparisonList = []
				count = count + 1
			else:
				rankComparisonList.append(tempRankList)
				tempRankList = [] # resets the list
			while (j < num4):
				if (i != j): 
					# If the two genes are connected, hasEdge = 1
					if (G.has_edge(geneIDList[i],geneIDList[j])): # this will need to change as will the below check case
						hasEdge = 1
					else:
						hasEdge = 0
					# temp_connection is the wij rj[n-1] / degi part, for this gene 
					# The if statement is done so that should a gene have no connections
					# for any reason, the program won't break. The results however may
					# not be as expected.
					#if (G.degree(geneIDList[i]) == 0):
					#	temp_connection = (hasEdge * (rankingValueList[i-1])) / 1 
					#else:
					#	temp_connection = (hasEdge * (rankingValueList[i-1])) / (G.degree(geneIDList[i])) # need to get degree etc from adjaceny matrix if possible
					temp_connection = (hasEdge * (rankingValueList[i-1])) / (G.degree(geneIDList[i]))
					# sumOfConnection is the current sum of the above part, for i iterations
					sumOfConnectionList[j] = sumOfConnectionList[j] + temp_connection
					connectionValue = (1-d)*normExprDataList[j]
					rankingValueList[j] = connectionValue + d*(sumOfConnectionList[j]) # this forms the (1-d)exj part   # this takes the normalised ex value
					tempRankList.append(rankingValueList)
					
					if (count <= 1): 
						converged = False
					elif (bool(set(rankComparisonList[k]) & set(rankComparisonList[k-1]))):
						converged = True
					else:
						converged = False
					if (i == 1):
						print ("connectionValue:" , connectionValue)
						print ("SumValue:" , sumOfConnectionList[j])
						print ("num of itterations:" , numOfItterations)
				j = j + 1
		i = i + 1
		j = 0
		k = k + 1
	#print (gene_ex_r) 	# gene, expression value, normalised ex value, rank, sumOf, etc...
	print ("The algorithm converged:" , converged , "after" , numOfItterations , "itterations")
	return rankingValueList # need to return multiple items or print all to files to be read later.

# Creates a list for each gene and associated ranking value
# and adds this to a list containing all genes with their rankings.
# This list is in the form of [[Gene,RankingValue],[Gene,RankingValue]...]
# each gene is then printed with it's rank and ranking value. The list
# of all genes and ranking values is returned for later use.	
def sortAndPrintRanking(geneIDList, rankingValueList, geneNameList):
	geneAndRankList = []
	rankedList = []
	num = len(geneIDList)
	i = 0
	while (i < num):
		geneAndRankList.append(geneIDList[i])
		geneAndRankList.append(rankingValueList[i])
		geneAndRankList.append(geneNameList[i])
		rankedList.append(geneAndRankList)
		geneAndRankList = [] # reset list to null
		i = i + 1
		
	def getKey(item):
		return item[1]
	rankedList = sorted(rankedList, key = getKey, reverse = True)
	
	num2 = len(rankedList)
	i = 0
	while (i < num2):
		print(rankedList[i][0] , "(" , rankedList[i][2] , ")" , "is ranked:" , i+1, "with a ranking value of:", rankedList[i][1])
		i = i + 1
	return rankedList

def testValidity(rankedList):
	y_true = []
	y_score = []
	num6 = 1000	
	i = 0
	while (i < num6):
		tempValue = (rankedList[i][1])*1000
		y_true.append(tempValue)
		i = i + 1
		
	num7 = len(rankedList)	
	i = (len(rankedList) - 1000)
	while (i < num7):
		tempValue = (rankedList[i][1])*1000
		y_score.append(tempValue)
		i = i + 1
	y_true = np.array(y_true)
	#print (y_true)
	y_score = np.array(y_score)
	#print("Original ROC area: {:0.3f}".format(roc_auc_score(y_true, y_score)))
	print (metrics.auc(y_true, y_score))
	#rocScore = roc_auc_score(y_true, y_score)
	#print (rocScore)

fullGeneDataList = readFile();
G = makeGraph(fullGeneDataList);
G = connectGraph(G, fullGeneDataList);
geneIDList,exprDataList,normExprDataList,geneNameList = getMultipleLists(fullGeneDataList);
rankingValueList = geneRank(geneIDList, exprDataList, normExprDataList, G, 0.5);
sortAndPrintRanking(geneIDList, rankingValueList, geneNameList);

#rank = geneRank(geneIDList[A,B,C,D], exprData[0.25,.75,1.25,1.75], nomrmExprData[1,1,1,1], sumOfConnection[], G, d ):