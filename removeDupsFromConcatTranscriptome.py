

# first use the scripts getLongestIsoformsIllumina.py and labelFastaHeaders.py to get labelled, deIsoformed files for each tissue
# Then concatenate these files, and blastn them against each other, using no evalue threshold, and -outfmt '6 qseqid sseqid qlen slen length pident evalue'
# The output of the blastn file is then fed into here

def percentFloat(proportion,total):
	percentage = (float(proportion) * 100)/float(total)
	return(percentage)


def Cluster(blast):
	hitLists = []
	f = open(blast)
	blast = f.readlines()
	for i in blast:
		q = i.split()[0]
		d = i.split()[1]
		qlen = i.split()[2]
		slen = i.split()[3]
		length = i.split()[4]
		pident = i.split()[5]
		evalue = i.split()[6]
		pairLengths = list(qlen+slen)
		minPairLength = min(pairLengths)
		maxPairLength = max(pairLengths)
		percentLengthOfMinSeqLength = percentFloat(length,minPairLength)

		# limit to searches that arent self-blasted, have at least 99% identity and where the length of match is at least 90% of the smallest seqLength
		if q != d:
			if pident > 99:
				if percentLengthOfMinSeqLength > 90:
					hitPair = [q,d]
					hitLists.append(hitPair)
	# go through each pair and either merge with a another group, if one or more sequences are shared, or add to the list of groups. Write finished groups to a list
	merged = []
	outFile = open('concatTranscriptomeClusters', 'w')
	for lists in [set(d) for d in hitLists]:
		for m in merged:
			if not m.isdisjoint(lists):
				m.update(lists)
				break
		else:
			merged.append(lists)
	for i in merged:
		for s in i:
			outFile.write(s.rstrip('\n') + '\t')
		outFile.write('\n')


path = input('give me your blast file \n')
blastFile = (path)

Cluster(blastFile)


