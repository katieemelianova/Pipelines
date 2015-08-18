import subprocess

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

def findOrthologs(d):
	subprocess.call(["blastx -query " + d + " -db /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/TAIR10_pep_20101214_updated.fasta -evalue 1e-50 -out " + d.rstrip('fastaORF') + "BLASTX -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)
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
	ATprotDict = fastaDict('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/TAIR10_pep_20101214_updated.fasta')
	out = open(d.rstrip('fastaORF') + 'ATorthologs' ,'w')
	for i in AThits:
		out.write('>' + i + '\n' + ATprotDict[i] + '\n')

def ungap(string):
		f = open(string)
		fasta = f.read()
		outfile = open(string, 'w')
		fasta = fasta.replace('rframe-', 'rframe_m_')
		fasta = fasta.replace('-', '')
		fasta = fasta.replace('!', 'N')
		fasta = fasta.replace('_R_', '')
		outfile.write(fasta)		 
		f.close()
		outfile.close()

def getReadingFrame(A):
	f1 =open(A)
	aln = f1.readlines()
	print(aln)
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

	outfile = open(A + 'CorrectFrames', 'w')
	for i in range(len(seqNames)):
		b = dists[distList[i]]
		if not b.startswith('AT'):
			outfile.write('>' + b + '\n' + alnDict[b] + '\n')				
	outfile.close()

	ungap((A + 'CorrectFrames'))


path = input('give me a fasta file \n')
fasta = (path)
fastaSeqs = fastaDict(fasta)

outfile = open('seqTranslation.fasta', 'a')
for seq in fastaSeqs:
	tmpFile = open('fastaSeq', 'w')
	tmpFile.write('>' + seq + '\n' + fastaSeqs[seq] + '\n')
	tmpFile.close()
	findOrthologs('fastaSeq')
	subprocess.call(["python dna2pep-1.1/dna2pep.py -a -r all --outformat fasta fastaSeq > fastaSeqProt"], shell=True)
	subprocess.call(["cat fastaSeqProt fastaSeqATorthologs > fastaSeqAt"], shell=True)
	getReadingFrame('fastaSeqAt')
	f = open('fastaSeqAtCorrectFrames')
	out = f.read()
	outfile.write(out)







