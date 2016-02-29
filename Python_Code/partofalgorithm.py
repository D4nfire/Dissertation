connectedList = []
while (j < len(geneList)):	
	if (i != j):
		#print("This part works") # proof that the idea works
		#print(geneList[i])
		#print(geneList[j])
		if(bool(set(geneList[i]) & set(geneList[j]))):
			if(geneList[i] not in connectedList):
				connectedList.append(geneList[i])
				G.add_edge(geneList[i-1],geneList[j-1])
			elif(geneList[j] not in connectedList): 
				connectedList.append(geneList[j])
				G.add_edge(geneList[i-1],geneList[j-1])
						
			G.add_edge(geneList[i-1],geneList[j-1])
			print(geneList[i-1] + ", " + geneList[j-1])
			
			#print("success")  proof that this part works
			i = i + 2
			j = j + 2
	i = i + 2
	j = j + 2
	
	

while (j < len(geneList)):	
	if (i != j):
		#print("This part works") # proof that the idea works
		#print(geneList[i])
		#print(geneList[j])
		if(bool(set(geneList[i]) & set(geneList[j]))):			
			G.add_edge(geneList[i-1],geneList[j-1])
			print(geneList[i-1] + ", " + geneList[j-1])
			
			#print("success")  proof that this part works
			i = i + 2
			j = j + 2
	i = i + 2
	j = j + 2