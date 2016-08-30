# New pipeline without clustering and using OrthoFinder Clusters
# Clustered 4 species and have used the Fasta Header names to input into the pipeline_control file. e.g. Athaliana fasta headers are like 'AT1G55380.1' so AT is the identifier for Athaliana.
# I used protein sequences to cluster but will use transcript sequences for alignments etc.
# N/B. make sure there are no additional carriage returns or newlines in the pipeline_control file



import subprocess


################ translational alignment step

def codons(seq):
	# get longest sequence without a stop codon
	ORF_list = []
	codon_list = []
	for i in range(len(seq)):
		codon = [seq[j:j+3] for j in range(i, len(seq), 3)]
		
		current_seq	 = []
		for i in codon:
			if i != 'TAG' and i != 'TAA' and i != 'TGA':
				current_seq.append(i)
			else:
				joined_seq = ''.join(current_seq)
				ORF_list.append(joined_seq)
				del(current_seq[:])
				break			
	return(max(ORF_list, key = len))


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

def callORF(L):
	allORFs = {}
	for i in L:
		s = L[i]
		rcs = rc(s)
		stop = ('TAG','TAA','TGA')
		# check if the sequence is a fragment with no stop codons in it, if so, add it to the list to be returned
		if any(i not in s for i in stop) or any(i not in rcs for i in stop):
			name = (i)
			seq = (s)
			allORFs[name] = (seq)
		else:
			# if the sequence does have stop codons in it, use def codons to get the longest subsequence without a stop codon 
			FandR = []
			F = codons(s)
			R = codons(rc(s))
			FandR.append(F)
			FandR.append(R)
			name = (i)
			seq = max(FandR, key = len)
			allORFs[name] = (seq)
	return(allORFs)

def findOrthologs(d):
	# blasts against a TAIR database that you must have in a db somewhere (i.e. change paths as necessary)
	subprocess.call(["blastx -query " + d + " -db /Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_pep_20101214_updated.txt -max_target_seqs 1 -evalue 1e-40 -out " + d + "BLASTX -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)
	out = open(d + 'ATorthologs' ,'w')
	nucOut = open(d + 'ATnuclOrthologs' ,'w')
	f = open(d + 'BLASTX')
	topBlast = f.readlines()
	if topBlast:
		AThits = []
		for i in topBlast:
			beg = i.split()[0]
			at = i.split()[1]
			length = i.split()[2]
			pident = i.split()[3]
			evalue = i.split()[4]
			if int(length) > 200:
				if float(pident) > 40:
					AThits.append(at)
		AThits = set(AThits)



		ATprotDict = fastaDict('/Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_pep_20101214_updated.txt')
		ATnuclDict = fastaDict('/Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_seq_20101214_updated.txt')
		for i in AThits:
			nucOut.write('>' + i + 'OUTGROUP' + '\n' + ATnuclDict[i] + '\n')
			out.write('>' + i + 'OUTGROUP' + '\n' + ATprotDict[i] + '\n')

def rc(seq): 
	# Reverse complements a sequence
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-', 'W':'N', 'Y':'N', 'S':'N', 'K':'N', 'M':'N', 'R':'R', 'H':'H'} 
	useq = seq.upper()
	sort = str(type(useq))
	letters = list(useq) 
	completters = [basecomplement[base] for base in letters] 
	sort = str(type(completters))
	rc_list =  completters[::-1]
	sort = str(type(rc_list))
	rc_seq = "".join(rc_list)
	sort = str(type(rc_list))
	return rc_seq


def getReadingFrame(A, N):
	protOutfile = open(str(A) + 'correctFrames', 'w')
	nuclOutfile = open(str(N) + 'correctFramesNucl', 'w')
	subprocess.call(["clustalo --force -i " + A + " --distmat-out=" + A + ".dist --full"], shell=True)
	distList = []

	allFramesDict = fastaDict('myOrfsmafftProtAT')
	nucSeqsDict = fastaDict(N)

	f1 =open(A)
	aln = f1.readlines()
	seqNames = []
	for i in aln:
		if i.startswith('>'):
			i = i.split('_rframe')[0]
			if not 'OUTGROUP' in i:
				seqNames.append(i)
	seqNames = set(seqNames)

	f2 = open(A + ".dist")
	distMat = f2.readlines()
	dists = {}
	for i in distMat:
		item = i.split()
		name = (item[0])
		dist = (item[-1])
		distList.append(dist)
		dists[name] = (dist)
	for i in seqNames:
		seqDists = []
		i = i.lstrip('>')
		for d in dists:
			if d.startswith(i):
				seqDists.append(dists[d])
		for d in dists:
			if d.startswith(i):
				if 'OUTGROUP' not in d:
					if dists[d] == min(seqDists):
						seq = (allFramesDict[d])
						seq = seq.replace('-', '')
						protOutfile.write('>' + d + '\n' + seq + '\n')
						nuclOutfile.write('>' + (d.split('_rframe')[0]) + '\n' + nucSeqsDict[(d.split('_rframe')[0])] + '\n')

def checkOrientation(nuc, prot):
	codonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

	nucDict = fastaDict(nuc)
	protDict = fastaDict(prot)
	nucOut = open('clusters2pal2nalNT', 'w')
	protOut = open('clusters2pal2nalAA', 'w')
	for n in nucDict:
		nucName = (n)
		nucSeq = nucDict[n]
		for p in protDict:
			protName = (p)
			protSeq = protDict[p]
			if protName.startswith(nucName):
				if protName.endswith('rframe-1'):
					newSeq = (''.join(list(rc(nucSeq))[::-1]))
					newSeq = newSeq[::-1]
					nucOut.write('>' + nucName + '\n' + newSeq + '\n')
				elif protName.endswith('rframe-2'):
					newSeq = (''.join(list(rc(nucSeq))[::-1]))
					newSeq = newSeq[::-1]
					newSeq = newSeq[1:]
					nucOut.write('>' + nucName + '\n' + newSeq + '\n')
				elif protName.endswith('rframe-3'):
					newSeq = (''.join(list(rc(nucSeq))[::-1]))
					newSeq = newSeq[::-1]
					newSeq = newSeq[2:]
					nucOut.write('>' + nucName + '\n' + newSeq + '\n')
				elif protName.endswith('rframe1'):
					nucOut.write('>' + nucName + '\n' + nucSeq + '\n')
				elif protName.endswith('rframe2'):
					newSeq = nucSeq[1:]
					nucOut.write('>' + nucName + '\n' + newSeq + '\n')
				elif protName.endswith('rframe3'):
					newSeq = nucSeq[2:]
					nucOut.write('>' + nucName + '\n' + newSeq + '\n') 
				protOut.write('>' + protName + '\n' + protSeq + '\n')

def eliminateShortSeqs(seqFile):
	alnDict = fastaDict(seqFile)
	outfile = open(seqFile, 'w')
	for i in alnDict:
		if not 'OUTGROUP' in i:
			seq = alnDict[i]
			seq = seq.replace('-', '')
			if len(seq) > 300:
				outfile.write('>' + i + '\n' + alnDict[i] + '\n')

def translationalAlign(count):
	dictz = (fastaDict('cluster.fasta'))
	ORFdict = callORF(dictz)
	outfile = open('myOrfs','w')
	for o in ORFdict:
		outfile.write('>' + o + '\n' + ORFdict[o] + '\n')
	outfile.close()

	subprocess.call(["cat cluster.fastaATnuclOrthologs myOrfs > myOrfsAt"], shell=True)
	
	subprocess.call(["mafft --quiet --adjustdirection myOrfsAt > myOrfsAtMafft"], shell=True)
	mafftDict = fastaDict('myOrfsAtMafft')
	mafftOut = open('myOrfsMafft', 'w')
	for i in mafftDict:

		if not 'OUTGROUP' in i:
			ungappedSeq = (mafftDict[i]).replace('-', '')
			mafftOut.write('>' + i + '\n' + ungappedSeq + '\n')
	mafftOut.close()

	f = open('myOrfsMafft')
	fasta = f.read()
	seqCount = fasta.count('>')

	if seqCount > 1:	
		subprocess.call(["python dna2pep-1.1/dna2pep.py -a -r all --outformat fasta myOrfsMafft > myOrfsMafftProt"], shell=True)
		subprocess.call(["cat myOrfsMafftProt cluster.fastaATorthologs > myOrfsmafftProtAT"], shell=True)
		subprocess.call(["mafft --quiet --anysymbol myOrfsmafftProtAT > myOrfsmafftProtAT.aln"], shell=True)
		getReadingFrame('myOrfsmafftProtAT.aln', 'myOrfsMafft')
		checkOrientation('myOrfsMafftcorrectFramesNucl','myOrfsmafftProtAT.alncorrectFrames')
		subprocess.call(["mafft --anysymbol clusters2pal2nalAA > myOrfsmafftProtAT.alncorrectFrames.aln"], shell=True)
		subprocess.call(["perl /Users/katieemelianova/Desktop/Software/pal2nal.v14/pal2nal.pl -output fasta myOrfsmafftProtAT.alncorrectFrames.aln clusters2pal2nalNT > myOrfsPal2Nal" + str(count)], shell=True)
		subprocess.call(["java -jar /Users/katieemelianova/Desktop/Software/BMGE-1.12/BMGE.jar -i myOrfsPal2Nal" + str(count) + " -t DNA -of myOrfsPal2Nal" + str(count) + ".bmge"], shell=True)
		subprocess.call(["distmat -nucmethod Jukes-Cantor myOrfsPal2Nal" + str(count) + ".bmge -outfile myOrfsPal2Nal" + str(count) + ".distmat"], shell=True)
		

		

def translationalAlign2():
	dictz = (fastaDict('cluster.fasta'))
	ORFdict = callORF(dictz)
	outfile = open('myOrfs','w')
	for o in ORFdict:
		outfile.write('>' + o + '\n' + ORFdict[o] + '\n')
	outfile.close()

	subprocess.call(["cat cluster.fastaATnuclOrthologs myOrfs > myOrfsAt"], shell=True)
	
	subprocess.call(["mafft --quiet --adjustdirection myOrfsAt > myOrfsAtMafft"], shell=True)
	mafftDict = fastaDict('myOrfsAtMafft')
	mafftOut = open('myOrfsMafft', 'w')
	for i in mafftDict:

		if not 'OUTGROUP' in i:
			ungappedSeq = (mafftDict[i]).replace('-', '')
			mafftOut.write('>' + i + '\n' + ungappedSeq + '\n')
	mafftOut.close()
	#eliminateShortSeqs('myOrfsMafft')

	f = open('myOrfsMafft')
	fasta = f.read()
	seqCount = fasta.count('>')

	
	subprocess.call(["python dna2pep-1.1/dna2pep.py -a -r all --outformat fasta myOrfsMafft > myOrfsMafftProt"], shell=True)
	subprocess.call(["cat myOrfsMafftProt cluster.fastaATorthologs > myOrfsmafftProtAT"], shell=True)
	subprocess.call(["mafft --quiet --anysymbol myOrfsmafftProtAT > myOrfsmafftProtAT.aln"], shell=True)
	getReadingFrame('myOrfsmafftProtAT.aln', 'myOrfsMafft')
	checkOrientation('myOrfsMafftcorrectFramesNucl','myOrfsmafftProtAT.alncorrectFrames')
	subprocess.call(["mafft --anysymbol clusters2pal2nalAA > myOrfsmafftProtAT.alncorrectFrames.aln"], shell=True)
	subprocess.call(["perl /Users/katieemelianova/Desktop/Software/pal2nal.v14/pal2nal.pl -output fasta myOrfsmafftProtAT.alncorrectFrames.aln clusters2pal2nalNT > myOrfsPal2Nal"], shell=True)
	subprocess.call(["java -jar /Users/katieemelianova/Desktop/Software/BMGE-1.12/BMGE.jar -i myOrfsPal2Nal -t DNA -of myOrfsPal2Nal.bmge"], shell=True)
	subprocess.call(["distmat -nucmethod Jukes-Cantor myOrfsPal2Nal.bmge -outfile myOrfsPal2Nal.distmat"], shell=True)


################ end of translational alignment step


def getDistance(distFile, seqCount):
	distmat = open(distFile)
	lines = distmat.readlines()
	lines = (lines[8:])
	distsList = []
	for line in lines:
		line = line.rstrip('\n')
		dists = (line.split('\t'))
		dists = [d.lstrip() for d in dists if d]
		if seqCount > 2:
			dists = dists[1:-2]
			dists = [float(d) for d in dists if d]
			for d in (dists):
				distsList.append(d)
		elif seqCount == 2:
			if len(dists) == 3:
				dist = float(dists[1])
				distsList.append(dist)
	meanDist = sum(distsList)/len(distsList)
	maxDist = max(distsList)
	minDist = min(distsList)

	return meanDist, maxDist, minDist


def prepareDBs():
	# create a dict to hold all the dicts for all species
	dictionary_seq_dictionaries={}
	# Create an empty list to hold all species identifiers to refer to later
	identifiers_list=[]
	set_identifiers_list = []

	f = open('pipeline_control')
	testfile = f.readlines()
	# Read in contents of control file, creating a dict for each species, and adding these dicts to another dict
	for i in testfile:
		identifier=i.split()[0]
		identifiers_list.append(identifier)
		set_identifiers_list.append(identifier)
		path=i.split()[1]	
		dictionary_seq_dictionaries['{0}'.format(identifier)]=fastaDict(path)
		#filename.close()

	set_identifiers_list = [x.split('_')[0] for x in set_identifiers_list]
	return (identifiers_list, set_identifiers_list, dictionary_seq_dictionaries)

def getPairwiseExpressionCor(exp):
	comparedList = []
	pairwiseCorrelations = []
	if (len(exp)) > 1:
		for i in exp:
			for j in exp:
				pair1 = (i,j)
				pair2 = (j,i)
				if pair1 not in comparedList:
					if pair2 not in comparedList:
						if i != j:
							ExpressionCorrelation = []
							ExpressionCorrelation.append(i)
							ExpressionCorrelation.append(j)
							PairExpression = open('PairExpression', 'w')
							PairExpression.write('\t'.join((exp[i])[0:17])+ '\n')
							PairExpression.write('\t'.join((exp[j])[0:17])+ '\n')
							PairExpression.close()
							subprocess.call(["Rscript /Users/katieemelianova/Desktop/PipelineRunning/expressionCorrelation.R"], shell=True)
							PairExpressionCorrelation = open('/Users/katieemelianova/Desktop/PipelineRunning/PairExpressionCorrelation')
							ConPairExpressionCorrelation = PairExpressionCorrelation.read()
							ConPairExpressionCorrelation = ConPairExpressionCorrelation.rstrip('\n')
							ExpressionCorrelation.append(ConPairExpressionCorrelation)
							PairExpressionCorrelation.close()

							PairExpression = open('PairExpression', 'w')
							PairExpression.write('\t'.join((exp[i])[17:34])+ '\n')
							PairExpression.write('\t'.join((exp[j])[17:34])+ '\n')
							PairExpression.close()
							subprocess.call(["Rscript /Users/katieemelianova/Desktop/PipelineRunning/expressionCorrelation.R"], shell=True)	
							PairExpressionCorrelation = open('/Users/katieemelianova/Desktop/PipelineRunning/PairExpressionCorrelation')
							PlePairExpressionCorrelation = PairExpressionCorrelation.read()
							PlePairExpressionCorrelation = PlePairExpressionCorrelation.rstrip('\n')
							ExpressionCorrelation.append(PlePairExpressionCorrelation)
							PairExpressionCorrelation.close()

							PairExpression = open('PairExpression', 'w')
							PairExpression.write('\t'.join((exp[i])[0:17])+ '\n')
							PairExpression.write('\t'.join((exp[j])[17:34])+ '\n')
							PairExpression.close()
							subprocess.call(["Rscript /Users/katieemelianova/Desktop/PipelineRunning/expressionCorrelation.R"], shell=True)	
							PairExpressionCorrelation = open('/Users/katieemelianova/Desktop/PipelineRunning/PairExpressionCorrelation')
							ConPlePairExpressionCorrelation = PairExpressionCorrelation.read()
							ConPlePairExpressionCorrelation = ConPlePairExpressionCorrelation.rstrip('\n')
							ExpressionCorrelation.append(ConPlePairExpressionCorrelation)
							PairExpressionCorrelation.close()	
							pairwiseCorrelations.append(ExpressionCorrelation)


						comparedList.append(pair1)
						comparedList.append(pair2)
	return pairwiseCorrelations	

def RunPipeline(identifiers_list, set_identifiers_list, dictionary_seq_dictionaries):
	counter = 0
	
	# put expression values across all Becon sequences into a dict
	expression = open('TMMNormalizedCountsNoCR2NoPP3')
	expressions = expression.readlines()
	expressionDict = {}
	for e in expressions:
		transcriptName = e.split()[0]
		expressionValues = e.split()[1:36]
		expressionDict[transcriptName] = expressionValues


	f = open('OrthologousGroupsCHS.txt')
	outfile = open('clusterOutFile', 'w')
	clusters = f.readlines()
	for cluster in clusters:
		outfile.write('// \n\n')
		counter+=1
		# define dicts and lists to store cluster specific values to, and open file to write sequences of cluster to a FASTA file
		copyDict = {}
		seq_count = []
		filtered_seq_count = 0
		seq_list = []
		file = open('cluster.fasta', 'w')

		outfile.write('Copy Numbers: ')
		#Count taxon specific copy numbers and store counts in a list
		for x in set_identifiers_list:
			count = cluster.count(x)
			copyDict[x] = count
			seq_count.append(count)
			outfile.write(str(count) + '\t')
		outfile.write('\n\n')
		seq = cluster.split('\t')

		# get pairwise expression values
		clusterExpressionDict = {}
		for s in seq:
			if s.startswith('Becon'):
				if s in expressionDict:
					seqExpression = expressionDict[s]
					clusterExpressionDict[s] = seqExpression
		
		pairwiseExpressionCorrelation = getPairwiseExpressionCor(clusterExpressionDict)
		if (len(pairwiseExpressionCorrelation)) > 0:
			outfile.write('Start pairwise correlations (gene pair, ConDifference, PleDifference, conVpleDifference) \n')
			for p in pairwiseExpressionCorrelation:
				outfile.write('\t'.join(p) + '\n')
			outfile.write('End pairwise correlations \n\n')

		outfile.write('Sequence Names: ')
		for i in seq:
			seq_list.append(i)
			outfile.write(i + '\t')
		outfile.write('\n')
		#write cluster sequences to FASTA file
		outfile.write('Start Sequence Fasta \n')
		for i in identifiers_list:
			for x in seq:
				if i != '\n':
					if x != '\n':
						if x.startswith(i):
							x = x.rstrip('.p\n')

							if len((dictionary_seq_dictionaries[i])[x]) > 300:
								if x.startswith('Becon'):
								
									file.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')
									filtered_seq_count+=1
							outfile.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')
		outfile.write('End Sequence Fasta \n\n')
		file.close()
		if filtered_seq_count > 1:
			if filtered_seq_count < 30:


				findOrthologs('cluster.fasta')
				f = open('cluster.fastaATnuclOrthologs')
				orthologs = f.readlines()
				if orthologs:
					f = open('ATH_GO_GOSLIM.txt')
					GO = {}
					goFile = f.readlines()
					for i in goFile:
						AT = i.split('\t')[0]
						ATGO = i.split('\t')[4]
						ATGO = ATGO.lstrip()
						ATGO = ATGO.rstrip()
						GO[AT] = ATGO
					
					atOrthologList = []
					outfile.write('Athaliana Orthologs: ')
					for i in orthologs:
						if i.startswith('>'):
							atOrtholog = i.split('.')[0]
							atOrtholog = atOrtholog.lstrip('>')
							atOrtholog = atOrtholog.rstrip()
							atOrthologList.append(atOrtholog)
							outfile.write(str(atOrtholog) + '\t')
					outfile.write('\n')#
					outfile.write('GO terms: ')
					for a in atOrthologList:
						try:
							outfile.write(GO[a] + '\t')
						except KeyError:
							continue
					outfile.write('\n')
					
					#eliminateShortSeqs('cluster.fasta')
					translationalAlign2()
					
					try:
						meanDistance, maxDistance, minDistance = getDistance('myOrfsPal2Nal.distmat', filtered_seq_count)
						outfile.write('\nStart Pairwise Distance Stats (mean, max, min) \n')
						outfile.write(str(meanDistance) + '\t' + str(maxDistance) + '\t' + str(minDistance) + '\n')
						outfile.write('End Pairwise Distance Stats \n')
					except ZeroDivisionError:
						outfile.write('\nStart Pairwise Distance Stats (mean, max, min) \n')
						outfile.write('Pairwise distance unavailable \n')
						outfile.write('End Pairwise Distance Stats \n')


	outfile.write('\n')
	


idList, setidList, dictSeqDict = prepareDBs()
RunPipeline(idList, setidList, dictSeqDict)





