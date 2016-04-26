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

# This file will give the ranking of the KO gene for all values of d for each of
# the 40 test/validation files. The file names are hard-cooded as this is purely
# for testing not use by users. The GOTermIndex value should be changed according
# to the desired network for different tests.

# Read in and parse data from a given file. This file must be in a certain format 
# where geneID must be the first element of each line, expr data must be the 5th element,
# gene name must be the 7th element and the GO id's element(s) can be determined by the 
# user. This information is put into a list called fullGeneDataList.
# genes with their GO 
def readFile(filePath):
	#Reading and parsing the file
	fullGeneDataList = []
	tempList = []
	updatedExpr = 0

	with open(filePath) as infile:
		next(infile)
		for line in infile:
			tempList = line.split("	")
			# if the line doesn't contain the gene name or expr value then skip that line.
			if ((tempList[6] == "") or (tempList[4] == "")):
				next(infile)
			# otherwise parse the data to fullGeneDataList
			elif (tempList[6] not in fullGeneDataList):
				fullGeneDataList.append(tempList[6]) # Gene ID
				#fullGeneDataList.append(abs(float(tempList[4])))
				fullGeneDataList.append(abs(float(tempList[4])))
				fullGeneDataList.append(1)
			else: 
				num = len(fullGeneDataList)
			
				i = 0
				while (i < num):
					if (fullGeneDataList[i] == tempList[6]):
						updatedCount = fullGeneDataList[i+2] + 1
						fullGeneDataList[i+2] = updatedCount
						#updatedExpr = fullGeneDataList[i+1] + (abs(float(tempList[4])))
						updatedExpr = fullGeneDataList[i+1] + (abs(float(tempList[4])))
						fullGeneDataList[i+1] = updatedExpr
					i = i + 1
					
	num = len(fullGeneDataList)
	i = 0
	while (i < num):
		if (fullGeneDataList[i+2] > 0):
			updatedExpr  = fullGeneDataList[i+1] / fullGeneDataList[i+2]
			fullGeneDataList[i+1] = updatedExpr
		i = i + 3
		
	# This list has the form of [Gene ID, GO ID list, expr data, gene name, Gene ID...] 
	return fullGeneDataList

# Takes in three lists, geneIDList, rankingValueList and geneNameList, 
# combines the data into one list, rankedList which has the form of
# [geneID, ranking value, gene name, geneID...] This new list is then
# sorted based on the ranking value, so the first gene in the list will
# be the one with the highest ranking value.
def sortByRanking(fullGeneDataList):
	geneAndRankList = []
	rankedList = []
	num = len(fullGeneDataList)
	
	i = 0
	while (i < num):
		geneAndRankList.append(fullGeneDataList[i])
		geneAndRankList.append(fullGeneDataList[i+1])
		rankedList.append(geneAndRankList) 
		geneAndRankList = [] # reset list to null for the next genes data
		i = i + 3
	
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
def writeResultsToFile(rankedList, koGene, sumOfRanks, outputFile):
	num2 = len(rankedList)
	this_ranking = ""
	this_string = ""
	
	with open(outputFile, "a") as this_file:
	
		# Write out the ranking data for each gene
		i = 0
		while (i < num2):
			if (str(rankedList[i][0]) == koGene):
				this_ranking = [str(rankedList[i][0]), ", ", str(i+1), "\n"]
				this_file.writelines(this_ranking)
				sumOfRanks = sumOfRanks + i + 1
			i = i + 1
	
		#this_string = ["\nThe roc score is: " , str(rocScore), "\n\n\n\n"]
		#this_file.writelines(this_string)
		this_file.close()
	
	return sumOfRanks

# Write the header information to a specified file, then close the file
# for safety. This appends to file and does not override the files contents.
def writeHeaderToFile(outputFile):
	this_string = ""
	
	with open(outputFile, "a") as this_file:
		this_string = ["GeneName, rank based on expr value \n"]
		this_file.writelines(this_string)
	this_file.close()

def writeAverageRankingToFile(sumOfRanks,outputFile):
	sumOfRanks = sumOfRanks/40 
	this_string = ""
	
	with open(outputFile, "a") as this_file:
		this_string = ["Average rank:", str(sumOfRanks), "\n"]
		this_file.writelines(this_string)
	this_file.close()
	
def main(filePath, koGene, sumOfRanks, outputFile):
	fullGeneDataList = readFile(filePath);
	rankedList = sortByRanking(fullGeneDataList);
	sumOfRanks = writeResultsToFile(rankedList, koGene, sumOfRanks, outputFile);
	
	return sumOfRanks
	
writeHeaderToFile("Rankings.txt");	

# Hard-coded values for graph file and output file names as well
# as a list of all KO genes to be tested. Also hard-coded is the
# index for the GO ID's in the given input files. This is done to
# quickly run the ranking for all 40 genes for all values of d for 
# testing and validation purposes and is not designed for end users.
genePath = "C:\ThirdYear\Dissertation\Python_Code\GR_V1.9\DataFiles\\"
testGeneList = ["Abca1", "Btk", "Cav1", "Cav3", "Cftr", "Clcn1", "Cnr1", "Emd", "Epas1", "Esrra", "Gap43", "Gnmt", "Hdac1", "Hdac2", "Hsf4", "Hspa1a", "Il6", "Lhx1", "Lhx8", "Lmna", "Mbnl1", "Mst1r", "Myd88", "Nos3", "Phgdh", "Pmp22", "Ppara", "Pthlh", "Rab3a", "Rasgrf1", "Rbm15", "Runx2", "Scd1", "Slc26a4", "Srf", "Tcf7", "Tgm2", "Zc3h12a", "Zfp36", "Zfx"]
endOfPath = "_100nn.tsv"
sumOfRanks = 0

num = 40
# Runs main for each KO gene with hard-coded inputs
i = 0
while (i < num):
	geneAndAllRankingsList = []
	testCount = 0
	CountPerGene = 0
	sumOfRanks = main(genePath + testGeneList[i] + endOfPath, testGeneList[i], sumOfRanks, "Rankings.txt");
	i = i + 1

writeAverageRankingToFile(sumOfRanks, "Rankings.txt");