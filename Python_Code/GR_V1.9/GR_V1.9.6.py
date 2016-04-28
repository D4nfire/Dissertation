def readFile(filePath):
	#Reading and parsing the file
	temp_list = []
	rankingDataList = []
	topTenCountList = []
	topTwentyCountList = []
	count = 0
	temp_sum_value = 0
	temp_top_ten_count = 0
	temp_top_twenty_count = 0
	
	with open(filePath) as infile:
		next(infile)
		for line in infile:
			temp_list = line.split(',')
			# if the line is the header line then skip
			if (line.startswith( 'GeneName' ) is True):
				next(infile)
			# otherwise parse the data to fullGeneDataList
			elif (count == 0):
				i = 1
				while (i < 22):
					rankingDataList.append(temp_list[i])
					topTenCountList.append(0)
					topTwentyCountList.append(0)
					i = i + 1
				count =  count + 1
			else:
				i = 0
				while (i < 21):
					temp_sum_value = float(rankingDataList[i]) + float(temp_list[i+1])
					rankingDataList[i] = temp_sum_value
					if (float(temp_list[i+1]) <= 10):
						temp_top_ten_count = topTenCountList[i] + 1
						topTenCountList[i] = temp_top_ten_count
						
					if (float(temp_list[i+1]) <= 20):
						temp_top_twenty_count = topTwentyCountList[i] + 1
						topTwentyCountList[i] = temp_top_twenty_count
					i = i + 1
					
	print (rankingDataList)
	
	num = len(rankingDataList)
	i = 0
	while (i < num):
		temp_avg_value = float(rankingDataList[i]) / 38 # 38 for p-p, 40 for others
		rankingDataList[i] = temp_avg_value
		i = i + 1
		
	return rankingDataList, topTenCountList, topTwentyCountList
	
def writeToFile(outputFile, rankingDataList, topTenCountList, topTwentyCountList):
	with open(outputFile, "a") as this_file:
		num = len(rankingDataList)
		
		this_file.writelines("\n\nAverage: ")
		
		i = 0
		while (i < num):
			this_string = [str(rankingDataList[i]), " "]
			this_file.writelines(this_string)
			i = i + 1
		
		this_file.writelines("\ntop10 count: ")
		
		i = 0
		while (i < num):
			this_string = [str(topTenCountList[i]), " "]
			this_file.writelines(this_string)
			i = i + 1
		
		this_file.writelines("\ntop20 count: ")
		
		i = 0
		while (i < num):
			this_string = [str(topTwentyCountList[i]), " "]
			this_file.writelines(this_string)
			i = i + 1
		
		this_file.close()

rankingDataList = []
filepath = ""
outputFile = ""

filepath = input('Give the full file path: ')
outputFile = input('What do you want to call the output file? ')	

rankingDataList, topTenCountList, topTwentyCountList = readFile(filepath);
writeToFile(outputFile, rankingDataList, topTenCountList, topTwentyCountList);