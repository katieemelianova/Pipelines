# process file called 'CONtairCDNAe50' with outfmt = qseqid sseqid length qlen slen pident evalue 

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

def percentFloat(proportion,total):
	percentage = (float(proportion) * 100)/float(total)
	return(percentage)


begoniaDict = fastaDict('CON_unique_fasta')
outfile = open('CONuniqueCompleteish', 'a')

f = open('CONtairCDNAe50')
blast = f.readlines()
for b in blast:
	q = b.split()[0]
	d = b.split()[1]
	length = b.split()[2]
	qlen = b.split()[3]
	slen = b.split()[4]
	percQueryDatabase = percentFloat(float(qlen), float(slen))
	percLengthDatabase = percentFloat(float(length), float(slen))
	if percQueryDatabase > 60:
		outfile.write('>' + q + '\n' + begoniaDict[q] + '\n')

