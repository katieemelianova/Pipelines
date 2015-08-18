import subprocess



def blastn(q,d,e,o):
	subprocess.call(["blastn -query " + q + " -db " + d + " -evalue " + e + " -out " + o + " -outfmt '6 qseqid sseqid length qlen slen pident evalue'"], shell=True)

def getGO(F):
	f = open(F)
	out = open('blast2GO', 'w')
	lines = f.readlines()
	for i in lines:
		at = i.split()[1]
		at = at.split('.')[0]
		try:
			out.write(at + '\t' + goDict[at] + '\n')
		except KeyError:
			continue
blastn('CONuniqueCompleteish', 'TAIR10_cdna_20101214_updated.txt', '50', 'blastOut')


goDict = {}
f = open('ATH_GO_GOSLIM.txt')
GO = f.readlines()
for g in GO:
	at = g.split('\t')[0]
	goTerm = g.split('\t')[4]
	goDict[at] = goTerm
	f.close()

getGO('blastOut')