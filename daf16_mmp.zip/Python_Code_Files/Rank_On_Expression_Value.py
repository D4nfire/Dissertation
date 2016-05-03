# web sources for some parts of the code. 		

# Sorting a list in python used in sortAndPrintRanking:
# http://pythoncentral.io/how-to-sort-a-list-tuple-or-object-with-sorted-in-python/

# This file will give the ranking of the KO gene based on the absolute
# or actual expression change values for each gene. The file names are 
# hard-cooded as this is purely for testing and not for users. The GOTermIndex 

# Read in and parse data from a given file. This file must be in a certain format 
# where geneID must be the first element of each line, expr data must be the 5th element,
# gene name must be the 7th element. 
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
				#fullGeneDataList.append((float(tempList[4]))) # for actual expression change
				fullGeneDataList.append(abs(float(tempList[4]))) # for absolute expression change
				fullGeneDataList.append(1) # the number of times the gene appears so far
			# If the gene is already in the list then update the count showing
			# how many times the gene has been found and add chnage the 
			# expression change value to be the sum of this value for the multiple
			# instances found.
			else: 
				num = len(fullGeneDataList)
			
				i = 0
				while (i < num):
					if (fullGeneDataList[i] == tempList[6]):
						updatedCount = fullGeneDataList[i+2] + 1
						fullGeneDataList[i+2] = updatedCount
						#updatedExpr = fullGeneDataList[i+1] + ((float(tempList[4])))
						updatedExpr = fullGeneDataList[i+1] + (abs(float(tempList[4])))
						fullGeneDataList[i+1] = updatedExpr
					i = i + 1
	
	# update the expression value to be the average value for genes
	# which appear multiple times.
	num = len(fullGeneDataList)
	i = 0
	while (i < num):
		if (fullGeneDataList[i+2] > 0):
			updatedExpr  = fullGeneDataList[i+1] / fullGeneDataList[i+2]
			fullGeneDataList[i+1] = updatedExpr
		i = i + 3
		
	# This list has the form of [Gene ID, GO ID list, expr data, gene name, Gene ID...] 
	return fullGeneDataList

# Puts the gene names and expression values for each gene in a new list and
# then sorts this list based on expression value.
def sortByRanking(fullGeneDataList):
	geneAndRankList = []
	rankedList = []
	num = len(fullGeneDataList)
	
	# create rankedList
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
	
	# rankedList has the form [[gene name, expr value], [gene name, expr value]... etc]
	return rankedList

# Write the header information to a specified file, then close the file
# for safety. This appends to file and does not override the files contents.
def writeHeaderToFile(outputFile):
	this_string = ""
	
	with open(outputFile, "a") as this_file:
		this_string = ["GeneName, rank based on expr value \n"]
		this_file.writelines(this_string)
	this_file.close()	
	
# Takes in the ranking list, the KO gene, a count value
# and the desired output file name. Then writes to a file
# the KO gene name and the associated ranking.
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
	
		this_file.close()
	
	return sumOfRanks

# writes to a file the average rank of all 40 KO genes
def writeAverageRankingToFile(sumOfRanks,outputFile):
	sumOfRanks = sumOfRanks/40 
	this_string = ""
	
	with open(outputFile, "a") as this_file:
		this_string = ["Average rank:", str(sumOfRanks), "\n"]
		this_file.writelines(this_string)
	this_file.close()

# The main method which calls the other required methods 
def main(filePath, koGene, sumOfRanks, outputFile):
	fullGeneDataList = readFile(filePath);
	rankedList = sortByRanking(fullGeneDataList);
	sumOfRanks = writeResultsToFile(rankedList, koGene, sumOfRanks, outputFile);
	
	return sumOfRanks

# only want this to happen once	
writeHeaderToFile("Rankings.txt");	

# Hard-coded values for the output file names as well
# as a list of all KO genes to be tested. 
genePath = "C:\ThirdYear\Dissertation\daf16_mmp.zip\DataFiles\\"
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

# want this to happen once at the end
writeAverageRankingToFile(sumOfRanks, "Rankings.txt");