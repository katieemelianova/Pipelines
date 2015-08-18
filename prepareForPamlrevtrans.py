import os, sys, subprocess


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

def ungapAA():
	for i in dirs:
		if 'AA' in i:
			f = open(i)
			fasta = f.read()
			outfile = open(i + '_ungapped', 'w')
			fasta = fasta.replace('-', '')
			fasta = fasta.replace('!', 'N')
			outfile.write(fasta)
			f.close()
			outfile.close()

def ungap(string):
	dirs = os.listdir(path)
	for i in dirs:
		if string in i:
			f = open(i)
			fasta = f.read()
			outfile = open(i, 'w')
			fasta = fasta.replace('rframe-', 'rframe_m_')
			fasta = fasta.replace('-', '')
			fasta = fasta.replace('!', 'N')
			fasta = fasta.replace('_R_', '')
			outfile.write(fasta)		 
			f.close()
			outfile.close()
		
def findOrthologs(d):
	subprocess.call(["blastx -query " + d + " -db /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/TAIR10_pep_annotations -evalue 1e-50 -out " + d.rstrip('fastaORF') + "BLASTX -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)
	f = open(d.rstrip('fastaORF') + 'BLASTX')
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

	ATprotDict = fastaDict('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/TAIR10_pep_annotations')
	ATnuclDict = fastaDict('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/TAIR10_cdna_annotations')
	out = open(d.rstrip('fastaORF') + 'ATorthologs' ,'w')
	nucOut = open(d.rstrip('fastaORF') + 'ATnuclOrthologs' ,'w')
	for i in AThits:
		nucOut.write('>' + i + '\n' + ATnuclDict[i] + '\n')
		out.write('>' + i + '\n' + ATprotDict[i] + '\n')

def cleanup():
	subprocess.call(["rm *ungappedORFMafft *ungappedORF *NT.fasta_ungapped"], shell=True)


def rc(seq): 
	# Reverse complements a sequence
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-', 'W':'N', 'Y':'N', 'S':'N', 'K':'N', 'M':'N'} 
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


def reorderAlignment(n, p):
	nuclDict = fastaDict(n)
	protDict = fastaDict(p)
	protOutfile = open((p.rstrip('CorrectFrames.Aln') + '.protOrdered'), 'w')
	nucOutfile = open((n.rstrip('Mafft') + 'nucOrdered'), 'w')
	for i in protDict:
		protOutfile.write('>' + i + '\n' + protDict[i] + '\n')
		i = i.split('_rframe')[0]
		nucOutfile.write('>' + i + '\n' + nuclDict[i] + '\n')

def getReadingFrame(A, N):
	f1 =open(A)
	aln = f1.readlines()
	subprocess.call(["clustalo --force -i " + A + " --distmat-out=" + A + ".dist --full"], shell=True)
	f2 = open(A + ".dist")
	distMat = f2.readlines()
	#get all unique sequence names - number of hits to look for later
	seqNames = []
	for i in aln:
		if i.startswith('>'):
			i = i.split('_rframe')[0]
			seqNames.append(i)
	seqNames = set(seqNames)
	f1.close()
	#make the alignment into a dict of name and sequence
	alnDict = fastaDict(A)
	dists = {}
	distList = []
	for i in distMat:
		item = i.split()
		name = (item[0])
		dist = (item[-1])
		dists[dist] = (name)
		distList.append(dist)
	distList.sort()
	nuclDict = fastaDict(N)
	nucOutfile = open(N + 'RC', 'w')

	outfile = open((A.rstrip('protATaln')) + 'CorrectFrames', 'w')
	for i in range(len(seqNames)):
		b = dists[distList[i]]
		if not b.startswith('AT'):
			outfile.write('>' + b + '\n' + alnDict[b] + '\n')				
	outfile.close()

	ungap(((A.rstrip('protATaln')) + 'CorrectFrames'))




path = "./"
dirs = os.listdir(path)

#for i in dirs:
#	if 'WorkingSeqs' in i:
#		subprocess.call(["java -jar /Volumes/BACKUP/Bioinformatics/macse_v1.01b.jar -prog alignSequences -seq " + i], shell=True)







ungap('NT')
dirs = os.listdir(path)
for i in dirs:
	if 'NT' in i:
		dict = (fastaDict(i))
		ORFdict = callORF(dict)
		outfile = open(i + 'ORF','w')
		for o in ORFdict:
			outfile.write('>' + o + '\n' + ORFdict[o] + '\n')
		outfile.close()
dirs = os.listdir(path)
for d in dirs:
	if 'ORF' in d:
		findOrthologs(d)
dirs = os.listdir(path)
for i in dirs:
	if 'fastaORF' in i:
		subprocess.call(["cat  " + i.rstrip('fastaORF') + "ATnuclOrthologs " + i + " > " + i.rstrip('fastaORF') + "AT"], shell=True)
		subprocess.call(["mafft --adjustdirection --anysymbol " + i.rstrip('fastaORF') + "AT" + " > " + i.rstrip('fastaORF') + "Mafft"], shell=True)


dirs = os.listdir(path)
ungap('Mafft')

for i in dirs:
	if 'Mafft' in i:
		dictz = fastaDict(i)
		outfile = open(i, 'w')
		for i in dictz:
			if not i.startswith('AT'):
				outfile.write('>' + i + '\n' + dictz[i] + '\n')
		outfile.close()

dirs = os.listdir(path)
for i in dirs:
	if 'Mafft' in i:
		subprocess.call(["python dna2pep-1.1/dna2pep.py -a -r all --outformat fasta " + i + "  > " + (i.rstrip('Mafft')) + "prot"], shell=True)
		# NB always make sure that the AT ortholog is added at the END!!
		subprocess.call(["cat " + (i.rstrip('fastaMafft')) + "prot " + (i.rstrip('fastaMafft')) + "ATorthologs > " + i.rstrip('fastaMafft') + "protAT"], shell=True)
		subprocess.call(["mafft --anysymbol " + i.rstrip('Mafft') + "protAT > " + i.rstrip('Mafft') + "protATaln"], shell=True)
		protFilename = (i.rstrip('Mafft') + "protATaln")
		getReadingFrame(protFilename, i)
		filename = (i.rstrip('Mafft') + 'CorrectFrames')
		subprocess.call(['mafft --anysymbol ' + filename + ' > ' + filename + '.Aln'], shell=True)
		reorderAlignment(i, (filename + '.Aln'))
		protFile = (i.rstrip('Mafft') + 'protOrdered')
		nucFile = (i.rstrip('Mafft') + 'nucOrdered')



		subprocess.call(['perl /Volumes/BACKUP/Bioinformatics/pal2nal.v14/pal2nal.pl ' +  protFile + ' ' + nucFile + ' > ' + i.rstrip('Mafft') + 'Pal2Nal' ], shell=True)


# this now works with all cloned CHS. Next, you need to test how PAML/codeml responds to some non-overlapping sequences. For this, get 454 complete sequences and cut bits out. Run PAML on both complete and fragmented and see difference. Try this with a conserved and divergent gene family.

