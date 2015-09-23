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
	subprocess.call(["blastx -query " + d + " -db /Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_pep_20101214_updated.txt -evalue 1e-50 -out " + d + "BLASTX -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)
	f = open(d + 'BLASTX')
	blast = f.readlines()
	AThits = []
	for i in blast:
		beg = i.split()[0]
		at = i.split()[1]
		length = i.split()[2]
		pident = i.split()[3]
		evalue = i.split()[4]
		if int(length) > 250:
			if float(pident) > 50:
				AThits.append(at)
	AThits = set(AThits)


	ATprotDict = fastaDict('/Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_pep_20101214_updated.txt')
	ATnuclDict = fastaDict('/Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_seq_20101214_updated.txt')
	out = open(d.rstrip('ORF') + 'ATorthologs' ,'w')
	nucOut = open(d.rstrip('ORF') + 'ATnuclOrthologs' ,'w')
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


def eliminateShortSeqs():
	alnDict = fastaDict('myOrfsMafft')
	outfile = open('myOrfsMafft', 'w')
	for i in alnDict:
		seq = alnDict[i]
		seq = seq.replace('-', '')
		if len(seq) > 300:
			outfile.write('>' + i + '\n' + alnDict[i] + '\n')


def countUngapped():
	try:
		alnDict = fastaDict('myOrfsPal2Nal')
		listDict = list(alnDict.values())
		seqLength = (len(listDict[0]))
		unGappedColumns = 0
		for i in range(seqLength):
			columnContents = []
			for s in alnDict:
				columnContents.append(alnDict[s][i])
			if '-' not in columnContents:
				unGappedColumns+=1
		return unGappedColumns
	except IndexError:
		print(alnDict)
		print('just printed broken ungapping def')



def prepareForPaml(count):
	dictz = (fastaDict('cluster.fasta'))
	ORFdict = callORF(dictz)
	outfile = open('myOrfs','w')
	for o in ORFdict:
		outfile.write('>' + o + '\n' + ORFdict[o] + '\n')
	outfile.close()

	findOrthologs('myOrfs')

	subprocess.call(["cat myOrfsATnuclOrthologs myOrfs > myOrfsAt"], shell=True)
	#56$%^&*()(*&^%$£$%^&*(*&^%$£))
	# SHOULD I SWITCH ROUND THE BLASTX FIND ORTHOLOGS AND THE FIND ORFS? WOULD BE MORE EFFECTIVE AT FINDING ORTHOLOGS THIS WAY
	
	subprocess.call(["mafft --quiet --adjustdirection myOrfsAt > myOrfsAtMafft"], shell=True)
	mafftDict = fastaDict('myOrfsAtMafft')
	mafftOut = open('myOrfsMafft', 'w')
	for i in mafftDict:

		if not 'OUTGROUP' in i:
			ungappedSeq = (mafftDict[i]).replace('-', '')
			mafftOut.write('>' + i + '\n' + ungappedSeq + '\n')
	mafftOut.close()
	eliminateShortSeqs()

	f = open('myOrfsMafft')
	fasta = f.read()
	seqCount = fasta.count('>')
	if seqCount > 1:	
		#print('CHECK ONE')
		subprocess.call(["python dna2pep-1.1/dna2pep.py -a -r all --outformat fasta myOrfsMafft > myOrfsMafftProt"], shell=True)
		#print('CHECK TWO')
		subprocess.call(["cat myOrfsMafftProt myOrfsATorthologs > myOrfsmafftProtAT"], shell=True)
		#print('CHECK THREE')
		subprocess.call(["mafft --quiet --anysymbol myOrfsmafftProtAT > myOrfsmafftProtAT.aln"], shell=True)
		#print('CHECK FOUR')
		getReadingFrame('myOrfsmafftProtAT.aln', 'myOrfsMafft')
		# AFTER THIS YOU CAN TAKE OUT THE AT ORTHOLOGS - MARK THEM FIRST TO BE REMOVED AFTER THIS
		#print('CHECK FIVE')
		

		checkOrientation('myOrfsMafftcorrectFramesNucl','myOrfsmafftProtAT.alncorrectFrames')
		#print('CHECK SIX')
		subprocess.call(["mafft --anysymbol clusters2pal2nalAA > myOrfsmafftProtAT.alncorrectFrames.aln"], shell=True)
		#print('CHECK SEVEN')
		subprocess.call(["perl /Users/katieemelianova/Desktop/Software/pal2nal.v14/pal2nal.pl -nogap -output fasta myOrfsmafftProtAT.alncorrectFrames.aln clusters2pal2nalNT > myOrfsPal2Nal"], shell=True)
		ungappedCols = countUngapped()
		if ungappedCols > 200:
			subprocess.call(["mv myOrfsPal2Nal myOrfsPal2Nal" + str(count)], shell=True)
		else:
			subprocess.call(["rm myOrfsPal2Nal"], shell=True)
		
################ end of translational alignment step














def blastn(file):
	# makes a BLAST directory and uses discontinuous megablast to BLAST concatenated file of all species against itself [NEEDS TO BE GENERALIZED]
	subprocess.call(["makeblastdb -in alltranscriptomes2blastn -dbtype nucl"], shell=True)
	subprocess.call(["blastn -num_threads 4 -task dc-megablast -query alltranscriptomes2blastn -db alltranscriptomes2blastn -out alltranscriptomes.blastnout -evalue 1e-40 -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)

def blastx(nuc, prot):
	# makes a BLAST directory and uses discontinuous megablast to BLAST concatenated file of all species against itself [NEEDS TO BE GENERALIZED]
	subprocess.call(["makeblastdb -in alltranscriptomes2blastx -dbtype prot"], shell=True)
	subprocess.call(["blastx -num_threads 4 -query alltranscriptomes2blastn -db alltranscriptomes2blastx -out alltranscriptomes.blastxout -evalue 1e-40 -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)
	subprocess.call(["rm alltranscriptomes2blast*"], shell=True)

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


def dictionary_seq(identifier, filename):
	# Returns a dict of a multi FASTA file with sequence header as keys and sequence as values.
	dictionary = {}
	for line in filename:
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
		dictionary[identifier + C_name] = C_seq
	return dictionary

def rc(seq): 
	
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-', '!':'!'} 
	useq = seq.upper()
	sort = str(type(useq))
	letters = list(useq) 
	completters = [basecomplement[base] for base in letters] 
	sort = str(type(completters))
	rc_list =  completters[::-1]
	sort = str(type(rc_list))
	rc_seq = "".join(rc_list)
	sort = str(type(rc_list))
	return(rc_seq)

def codons(seq):
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

def callORF(L):
	allORFs = {}
	for i in L:
		s = L[i]
		rcs = rc(s)
		stop = ('TAG','TAA','TGA')
		# if no stop codons in the sequence, print the sequence
		if any(i not in s for i in stop) or any(i not in rcs for i in stop):
			#print('>' + i + '\n' + s + '\n')
			name = (i)
			seq = (s)
			allORFs[name] = (seq)
		else:
			FandR = []
			F = codons(s)
			R = codons(rc(s))
			FandR.append(F)
			FandR.append(R)
			#print('>' + i + '\n' + max(FandR, key = len))
			name = (i)
			seq = max(FandR, key = len)
			allORFs[name] = (seq)
	return(allORFs)



























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
		identifier = identifier.upper()
		identifiers_list.append(identifier)
		set_identifiers_list.append(identifier)
		#identifier = identifier.split('_')[0]
		path=i.split()[1]	
		filename = open(path)
		dictionary_seq_dictionaries['{0}'.format(identifier)]=dictionary_seq(identifier, filename)
		filename.close()

	set_identifiers_list = [x.split('_')[0] for x in set_identifiers_list]
	set_identifiers_list = list(set(set_identifiers_list))


	# write contents of dict of dicts to a multifasta file, filtering by sequence length
	all_transN = open('alltranscriptomes2blastn', 'w')
	for i in identifiers_list:
		if '_N' in i:
			for x in (dictionary_seq_dictionaries[i]):
				if int(len((dictionary_seq_dictionaries[i])[x])) > 100:
					all_transN.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')

	# write contents of dict of dicts to a multifasta file, filtering by sequence length


	all_transP = open('alltranscriptomes2blastx', 'w')
	for i in identifiers_list:
		if '_P' in i:
			for x in (dictionary_seq_dictionaries[i]):
				if int(len((dictionary_seq_dictionaries[i])[x])) > 100:
					all_transP.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')
	return (set_identifiers_list, dictionary_seq_dictionaries)



def blastSearch():
	blastn('alltranscriptomes2blastn')
	blastx('alltranscriptomes2blastn', 'alltranscriptomes2blastx')
	subprocess.call(["cat alltranscriptomes.blastxout alltranscriptomes.blastnout > alltranscriptomes.blastout"], shell=True)
				
				
def Cluster():
	hitLists = []
	f = open('alltranscriptomes.blastnout')
	blast = f.readlines()
	for i in blast:
		q = i.split()[0]
		d = i.split()[1]
		l = i.split()[2]
		p = i.split()[3]
		e = i.split()[4]
		# limit to searches that arent self-blasted, have a length of > 150 and have > 40% identity. Add all match pairs with this description to a list of lists
		if q != d:
			if int(l) > 250:
				if float(p) > 60:
					hitPair = [q,d]
					hitLists.append(hitPair)
	# go through each pair and either merge with a another group, if one or more sequences are shared, or add to the list of groups. Write finished groups to a list
	merged = []
	outFile = open('clusters_L150P60', 'w')
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

def test():
	counter = 0
	f = open('clusters_L150P60')
	outfile = open('clusterTestOutFile', 'w')
	clusters = f.readlines()
	for cluster in clusters:
		counter+=1
		# define dicts and lists to store cluster specific values to, and open file to write sequences of cluster to a FASTA file
		copyDict = {}
		seq_count = []
		seq_list = []
		file = open('cluster.fasta', 'w')
		
		# Count taxon specific copy numbers and store counts in a list
		for x in set_identifiers_list:
			count = cluster.count(x)
			copyDict[x] = count
			seq_count.append(count)
			#print(seq_count)
		if sum(seq_count) > 1:
			# Store all sequence names in a list to be referred to later
			seq = cluster.split('\t')
			for i in seq:
				seq_list.append(i)
			# write cluster sequences to FASTA file
			for i in identifiers_list:
				if '_N' in i:
					for x in seq:
						if x.startswith(i.rstrip('_N')):
							file.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')
		file.close()
		prepareForPaml(counter)
		# I THINK I SHOULD PUT THE 'IF THERES AN AT IN THE CLUTSR' BEFORE CALLING PREPARE FOR PAML?
		#**************************&&^^&*(*&*(*&*(*&&&*(*&*(*&*())))))
		file.close()

		# write copy numbers and sequence names to file (only if there is an Athaliana ortholog associated with the family)
		if any(item.startswith('AT') for item in seq_list):
			copyNumbers = []
			clust = cluster.split('\t')
			for i in copyDict:
				copy = str(copyDict[i])
				copyNumbers.append(copy)
			copies = copyNumbers + clust
			copies = '\t'.join(copies)
			outfile.write(copies)


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
	identifier = identifier.upper()
	identifiers_list.append(identifier)
	set_identifiers_list.append(identifier)
	#identifier = identifier.split('_')[0]
	path=i.split()[1]	
	filename = open(path)
	dictionary_seq_dictionaries['{0}'.format(identifier)]=dictionary_seq(identifier, filename)
	filename.close()

set_identifiers_list = [x.split('_')[0] for x in set_identifiers_list]
set_identifiers_list = list(set(set_identifiers_list))


# write contents of dict of dicts to a multifasta file, filtering by sequence length
all_transN = open('alltranscriptomes2blastn', 'w')
for i in identifiers_list:
	if '_N' in i:
		for x in (dictionary_seq_dictionaries[i]):
			if int(len((dictionary_seq_dictionaries[i])[x])) > 300:
				all_transN.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')











#prepareDBs()
#blastSearch()
#Cluster()

test()



#Tuesday evening - changed the nuc seq checking frame def to work, now have switched everything into this script - doesnt seem to be working but going home now - RUN and TEST toorrow






	
		




