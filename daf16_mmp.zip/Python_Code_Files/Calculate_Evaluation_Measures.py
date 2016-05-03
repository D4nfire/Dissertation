# Reads in and parses data from the output file containing the
# ranking for all 40 KO genes.
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
			# otherwise, initialise the lists for the 21 values of d evaluated
			elif (count == 0):
				i = 1
				while (i < 22):
					rankingDataList.append(temp_list[i]) # initial sum of ranks
					topTenCountList.append(0) # initial count
					topTwentyCountList.append(0) # initial count
					i = i + 1
				count =  count + 1
			# otherwise parse the data 
			else:
				# for each value of d, update the measures.
				i = 0
				while (i < 21):
					temp_sum_value = float(rankingDataList[i]) + float(temp_list[i+1])
					rankingDataList[i] = temp_sum_value # update the sum of ranks
					# if the KO gene is in the top10 genes, add 1 to topTenCountList[i]
					if (float(temp_list[i+1]) <= 10):
						temp_top_ten_count = topTenCountList[i] + 1
						topTenCountList[i] = temp_top_ten_count 
					# if the KO gene is in the top10 genes, add 1 to topTwentyCountList[i]
					if (float(temp_list[i+1]) <= 20):
						temp_top_twenty_count = topTwentyCountList[i] + 1
						topTwentyCountList[i] = temp_top_twenty_count
					i = i + 1
	
	# calculate the average rank for each value of d	
	num = len(rankingDataList)
	i = 0
	while (i < num):
		temp_avg_value = float(rankingDataList[i]) / 40 # 38 for p-p network, 40 for all others
		rankingDataList[i] = temp_avg_value
		i = i + 1
	
	# return the lists containing the evaluation measures for each value of d
	return rankingDataList, topTenCountList, topTwentyCountList

# write the 3 lists to a file, usually the same as the input file.
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

# main methpd which asks for the input and output file names 
# and calls the required methods
def main():
	rankingDataList = []
	filepath = ""
	outputFile = ""

	filepath = input('Give the full file path: ')
	outputFile = input('What do you want to call the output file? ')	

	rankingDataList, topTenCountList, topTwentyCountList = readFile(filepath);
	writeToFile(outputFile, rankingDataList, topTenCountList, topTwentyCountList);

main();