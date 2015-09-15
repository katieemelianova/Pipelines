
import subprocess

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
		if int(length) > 100:
			if float(pident) > 70:
				AThits.append(at)
	AThits = set(AThits)


	ATprotDict = fastaDict('/Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_pep_20101214_updated.txt')
	ATnuclDict = fastaDict('/Users/katieemelianova/Desktop/Software/ncbi-blast-2.2.30+/db/TAIR10_seq_20101214_updated.txt')
	out = open(d.rstrip('ORF') + 'ATorthologs' ,'w')
	nucOut = open(d.rstrip('ORF') + 'ATnuclOrthologs' ,'w')
	for i in AThits:
		nucOut.write('>' + i + '\n' + ATnuclDict[i] + '\n')
		out.write('>' + i + '\n' + ATprotDict[i] + '\n')

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

	allFramesDict = fastaDict(A)
	nucSeqsDict = fastaDict(N)

	f1 =open(A)
	aln = f1.readlines()
	seqNames = []
	for i in aln:
		if i.startswith('>'):
			i = i.split('_rframe')[0]
			if not i.startswith('>AT'):
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


def prepareForPaml():
	dictz = (fastaDict(myFile))
	ORFdict = callORF(dictz)
	outfile = open('myOrfs','w')
	for o in ORFdict:
		outfile.write('>' + o + '\n' + ORFdict[o] + '\n')
	outfile.close()

	findOrthologs('myOrfs')

	subprocess.call(["cat myOrfsATnuclOrthologs myOrfs > myOrfsAt"], shell=True)
	subprocess.call(["mafft --adjustdirection myOrfsAt > myOrfsAtMafft"], shell=True)
	mafftDict = fastaDict('myOrfsAtMafft')
	mafftOut = open('myOrfsMafft', 'w')
	for i in mafftDict:
		if not i.startswith('AT'):
			ungappedSeq = (mafftDict[i]).replace('-', '')
			mafftOut.write('>' + i + '\n' + ungappedSeq + '\n')


	subprocess.call(["python dna2pep-1.1/dna2pep.py -a -r all --outformat fasta myOrfsMafft > myOrfsMafftProt"], shell=True)
	subprocess.call(["cat myOrfsMafftProt myOrfsATorthologs > myOrfsmafftProtAT"], shell=True)
	subprocess.call(["mafft myOrfsmafftProtAT > myOrfsmafftProtAT.aln"], shell=True)

	getReadingFrame('myOrfsmafftProtAT.aln', 'myOrfsMafft')

	subprocess.call(["mafft myOrfsmafftProtAT.alncorrectFrames > myOrfsmafftProtAT.alncorrectFrames.aln"], shell=True)

	subprocess.call(["perl /Users/katieemelianova/Desktop/Software/pal2nal.v14/pal2nal.pl -output fasta myOrfsmafftProtAT.alncorrectFrames.aln myOrfsMafftcorrectFramesNucl > myOrfsPal2Nal"], shell=True)





###########

PATH = input('Give me a fasta file \n')
myFile = (PATH)

prepareForPaml()


