

def fastaDict(i):
	f = open(i)
	L = {}

	for line in f:
		if line.startswith(">"):
			C_seq = ''
			C_split_line = line.split(' ')
			C_name = C_split_line[0]
			C_name = C_name.rstrip()
			C_name = C_name.lstrip('>')
		else:
			C_seq = C_seq + str(line)
			C_seq = C_seq.rstrip()
			C_seq = C_seq.upper()

		L[C_name] = C_seq
	return(L)

def countUngapped():
	alnDict = fastaDict('myOrfsPal2Nal2')
	listDict = list(alnDict.values())
	seqLength = (len(listDict[0]))
	unGappedColumns = 0
	for i in range(seqLength):
		columnContents = []
		for s in alnDict:
			columnContents.append(alnDict[s][i])
		if '-' not in columnContents:
			unGappedColumns+=1
	print(unGappedColumns)


countUngapped()
