
# This is done purely to get some expression value in
# to do the rest of the algorithm, the ex values here
# are randomly generated floats from -4 to +4. This
# is put into a list in the form of gene, ex, gene, ex etc
# This also initialises the ranking of each gene, making the
# list have the form of gene, ex, r, gene, ex, r etc...
i = 0
num = len(gene_ex_r)
num2 = 1
temp_ex = 0
sumOfEx = 0
while (i < num):
	temp_ex = random.uniform(-4.0, 4.0)
	gene_ex_r.insert(num2 , temp_ex)
	temp_ex = abs(temp_ex)
	sumOfEx = sumOfEx + temp_ex
	i = i + 1

temp_r = sumOfEx/(numpy.linalg.norm([sumOfEx], ord=1))
sumOfConnection = 0 # saves time later
i = 2
while (i < num):
	gene_ex_r.insert(num2 , temp_r)
	num2 = num2 + 1
	gene_ex_r.insert(num2, sumOfConnection)
	num2 = num2 + 2
#print (gene_ex_r) 	# gene, expression value, etc...

# Run num3 iterations, updating every genes rank
# per iteration
num3 = 100
num4 = len(gene_ex_r)
d = 0.5
connectionValue = 1 - d
i = 0 # getting string/int errors for parts later
j = 0
rank = 0
while (i < num3):
	print (i)
	while (j < num4):
		if (i != j):
			if (G.has_edge(gene_ex_r[i],gene_ex_r[j])):
				hasEdge = 1
			else:
				hasEdge = 0
			temp_connection = (hasEdge * (gene_ex_r[j+2])) / (G.degree(gene_ex_r[i]))
			sumOfConnection = (gene_ex_r[j+3]) + temp_connection
			gene_ex_r[i+3] = sumOfConnection
			rank = connectionValue + (d * gene_ex_r[j+3])
			gene_ex_r[i+2] = rank
		j = j + 3
	i = i + 3
	j = 0
print (gene_ex_r)



# rankJ = (1-d)*(gene_ex_r[J]) + (For all i (d*(1 or 0)*rankI*outDegree))
	
i = 0
j = 0
tempVal = 0
tempVal2 = 0
rank = 0
hasEdge = 0
while (i < num2):
	while (j < num2):
		if (i != j):
			if (G.has_edge(gene_ex_r[i],gene_ex_r[j])):
				hasEdge = 1
			else:
				hasEdge = 0
			tempVal = (tempVal + (hasEdge*(gene_ex_r[j+2])) / (G.degree(gene_ex_r[j])))
			j = j + 3
	tempVal = (tempVal * d)
	tempVal2 = ((1-d)*(gene_ex_r[i+1]))
	rank = tempVal + tempVal2
	gene_ex_r[i+2] = rank
	i = i + 3

print (gene_ex_r)