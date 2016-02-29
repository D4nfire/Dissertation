
systematicGeneList = []
averageExpressionList = []
tempExpressionList = []

#Open the yeast expression file and extract the gene name (systematic name)
#Then extract the expression data for the associated gene and take the 
#average of the data. Add to a list in format of Gene, Average Expression
#value, Gene, Average expression value, etc...
with open('C:\ThirdYear\Dissertation\Other helpful documents\Yeast_expression_data.txt') as infile:
	for line in infile:
		if (line.startswith('Y')):
			tempList = line.split("	")
			if (tempList[0] != "YORF"):
				systematicGeneName = tempList[0]
				systematicGeneList.append(systematicGeneName)
				#print(systematicGeneName) 		#Proof it works
				tempExpression = 0
				averageExpression = 0
				tempExpressionList = tempList[3:len(tempList)]
				tempExpressionList = map(float, tempExpressionList)
				tempExpression = sum(tempExpressionList)
				#print (tempExpression)		 #Proof this works
				numOfDataSamples = len(tempList) - 3
				averageExpression = tempExpression/numOfDataSamples
				#print (averageExpression)		#Proof this works
				#averageExpressionList.append(averageExpression)
				#print (averageExpressionList)		#proof this works
				systematicGeneList.append(averageExpression) # add the expression
				# data value to the list of genes, in format of gene, expression
				# value, gene, expression value etc...
				print (systematicGeneList)
				