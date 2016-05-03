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

# Read in and parse data from a given file. This file must be in a certain format 
# where geneID must be the first element of each line, expr data must be the 5th element,
# gene name must be the 7th element and the GO id's element(s) can be determined by the 
# user. This information is put into a list called fullGeneDataList.
def readFile(filePath, GOTermIndex):
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
				goList = (tempList[GOTermIndex].split("///")) #GO ID's
				fullGeneDataList.append(goList)
				fullGeneDataList.append(abs(float(tempList[4]))) # expr value
				fullGeneDataList.append(tempList[6]) # gene name
				goList = [] # resets the temporary list for the next line
	
	# This list has the form of [Gene ID, [GO ID list], expr data, gene name, Gene ID...] 
	return fullGeneDataList

# Read in and parse data from a given file. This file must be in a certain format 
# where geneID must be the first element of each line, expr data must be the 5th element,
# gene name must be the 7th element and the GO id's must be the 14th, 15th and 16th
# elements. This information is put into a list called fullGeneDataList.
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
				fullGeneDataList.append(tempList[0]) # add the Gene ID
				tempGoList1 = (tempList[13].split("///")) # GO.Function.ID
				tempGoList2 = (tempList[14].split("///")) # GO.Process.ID
				tempGoList3 = (tempList[15].split("///")) # GO.Component.ID
				goList = tempGoList1 + tempGoList2 + tempGoList3 # merge the three GO ID lists into one list
				fullGeneDataList.append(goList) # go ID list
				fullGeneDataList.append(abs(float(tempList[4]))) # expression change value (expr)
				fullGeneDataList.append(tempList[6]) # gene name
				goList = [] # resets the temporary list for the next line
	
	# This list has the form of [Gene ID, [GO ID list], expr data, gene name, Gene ID...] 
	return fullGeneDataList

# if a gene appears in the list multiple times we want to change the
# expr value to an average of the expr values for this gene and remove all
# later references to this gene. This is done by list manipulation.
def removeDuplicates(fullGeneDataList):
	duplicatedGeneList = []
	temp_gene_list = []
	num1 = len(fullGeneDataList)
	count = 0
	index = 0
	temp_expr = 0
	avgExprValue = 0
	
	# Get a list of duplicated genes, duplicatedGeneList
	i = 0
	while (i < num1):
		count = fullGeneDataList.count(fullGeneDataList[i+3])
		if ((count > 1) and (fullGeneDataList[i+3] not in duplicatedGeneList)):
			duplicatedGeneList.append(fullGeneDataList[i+3])
		i = i + 4
	
	# add two numbers for each gene in duplicatedGeneList, the firts number
	# will be the sum of expr values for this gene, the second will become
	# the number of times the gene is mentioned. These are set to 0 and 1 
	# initially
	# duplicatedGeneList will be in the form of [geneName, 0, 1, geneName, ...]
	num2 = len(duplicatedGeneList)*3
	i = 0
	while (i < num2):
		duplicatedGeneList.insert(i+1, 0)
		duplicatedGeneList.insert(i+2, 1)
		i = i + 3
	
	# This calculates and updates the sum of all expr values for the given duplicate genes as well
	# as the number of times the gene is found in the file.
	# duplicatedGeneList will be in the form of [geneName, sumOfExprForGene, numOfTimesInFile, ...]
	i = 0
	while (i < num1):
		# If the gene is a duplicated in the file 
		if (fullGeneDataList[i+3] in duplicatedGeneList):
			count = fullGeneDataList.count(fullGeneDataList[i+3]) # the number of times the gene is in the file
			index = duplicatedGeneList.index(fullGeneDataList[i+3]) + 1 # the index for the total expr value of the gene
			temp_expr = duplicatedGeneList[index] + fullGeneDataList[i+2] # update the total expr value for the gene
			duplicatedGeneList[index] =  temp_expr # update the total sum of expr values for the gene 
			duplicatedGeneList[index+1] = count # update the value showing how many times the gene is in the file
		i = i + 4
	
	# This changes the list element for the sum of all expr values for the gene
	# to be the average expr value for that gene.
	num3 = len(duplicatedGeneList)
	i = 0
	while (i < num3):
		avgExprValue = duplicatedGeneList[i+1] / duplicatedGeneList[i+2] # calculates the average expr value for the given gene
		duplicatedGeneList[i+1] = avgExprValue # updates this from sum of expr value to averageg expr value for the gene.
		i = i + 3
	
	# For duplicated genes update the first instance of the duplicate gene to have the
	# average expr value for the gene and remove/delete all other references to this gene
	i = 0
	while (i < num1):
		# If the gene is a duplicate and the first instance of the gene hasn't been updated 
		# to use the average expr value 
		if ((fullGeneDataList[i+3] in duplicatedGeneList) and (fullGeneDataList[i+3] not in temp_gene_list)):
			index = duplicatedGeneList.index(fullGeneDataList[i+3]) + 1 # get the average expr value for the gene
			fullGeneDataList[i+2] = duplicatedGeneList[index] # update the genes expr value to be the avg expr value
			temp_gene_list.append(fullGeneDataList[i+3]) # add the gene to a list to show it has been updated
			# Increase by 4 for the next gene
			i = i + 4
		
		# If the gene has had it's first instance updated to use the average expr value then
		# delete all duplicate gene information, so when the gene is found again, delete the
		# id, goList, exprValue and geneName of the duplicate. There is no need to update i
		# as the next gene will have moved down to i as the previous genes data was removed.
		elif (fullGeneDataList[i+3] in temp_gene_list):
			del fullGeneDataList[i] 
			del fullGeneDataList[i]
			del fullGeneDataList[i]
			del fullGeneDataList[i]
		
		# Otherwise increase i by 4 to get to the next gene
		else:
			i = i + 4
		
		# num1 must be updated to the new length of the list as removing items from
		# the list will shorten it.
		num1 = len(fullGeneDataList)
	
	return fullGeneDataList
	
# For each gene ID add a node to the graph G, return the graph for later us
# fullGeneDataList has the form of [Gene ID, GO ID list, expr data, gene name, Gene ID...] 
def makeGraph(fullGeneDataList):
	G = nx.Graph()
	
	i = 0
	while (i < len(fullGeneDataList)):
		G.add_node(fullGeneDataList[i])	# add a node for each gene
		i = i + 4
		
	return G

# For each gene, check the genes list of GO ID's against every other gene
# if two genes both share one or more GO ID's then add an edge between these
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
# Split the fullGeneDataList into four lists, geneIDList, exprDataList, normExprDataList
# and geneNameList as these are required for the main algorithm and later methods.
def getMultipleLists(fullGeneDataList):
	geneIDList = []
	exprDataList = []
	normExprDataList = []
	geneNameList = []
	sumOfEx = 0 # the sum of all expr values 
	
	i = 0
	j = 0
	while (i < len(fullGeneDataList)):
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
# as well as the Graph G and parameter d. The algorithm is run such that each gene is
# given an initial ranking score, the normalised expression value of that gene. The ranking
# score is updated each iteration. The algorithm runs for n iterations where n is the number 
# of genes. !The data in the three lists must be in the same gene order, such that 
# geneIDList[i], exprDataList[i] and normExprDataList[i] all belong to the same gene!
def geneRank(geneIDList, exprDataList, normExprDataList, G, d):
	sumOfConnectionList = []
	rankingValueList = []
	num = len(geneIDList)

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
				# sumOfConnection[j] is the current sum of the above part for gene j, after i iterations
				sumOfConnectionList[j] = sumOfConnectionList[j] + temp_connection
				connectionValue = (1-d)*exprDataList[j]
				rankValue = connectionValue + d*(sumOfConnectionList[j]) # The final ranking value
				rankingValueList[j] = rankValue # update the ranking score for gene j
				
			j = j + 1
		i = i + 1
		j = 0 # reset j so that every gene is compared again next iteration

	return rankingValueList # List of all genes ranking values in the order they were put in

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
	
	# rankedList has the form of [geneID, ranking value, gene name, geneID...]
	return rankedList

# rankedList has the form of [geneID, ranking value, gene name, geneID...]
# and koGene is the knocked out gene the experiment was based on. This method gives a
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
		if (rankedList[i][2] == koGene):
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
	
		this_string = ["Final ranking list and roc information for variable d = ", str(d), "\n\n"]
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
def runRanking(filepath, koGene, graphFile, outputFile, d, GOTermIndex):
	#fullGeneDataList = readFile(filepath, GOTermIndex);
	fullGeneDataList = readFile(filepath);
	fullGeneDataList = removeDuplicates(fullGeneDataList);
	G = makeGraph(fullGeneDataList);
	G = connectGraph(G, fullGeneDataList, graphFile);
	geneIDList,exprDataList,normExprDataList,geneNameList = getMultipleLists(fullGeneDataList);
	rankingValueList = geneRank(geneIDList, exprDataList, normExprDataList, G, d);
	rankedList = sortByRanking(geneIDList, rankingValueList, geneNameList);
	rocScore = testValidity(rankedList, koGene);
	writeResultsToFile(rankedList, rocScore, d, outputFile);

# Runs the file and asks for the file path information and the KO gene
# name. d is set to 0.65 but can be cahnged.
def main():
	filepath = input('Give the full file path: ')
	koGene = input('Give the name of the ko gene: ')
	graphFile = input('What do you want to call the graph file? include .xml ') 
	outputFile = input('What do you want to call the output file? ')
	#GOTermIndex = input('What is the index of the GO term?')
	
	runRanking(filepath, koGene, graphFile, outputFile, 0.65, "");
	#runRanking(filepath, koGene, graphFile, outputFile, 0.65, GOTermIndex);
	
main();