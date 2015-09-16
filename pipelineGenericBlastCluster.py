import subprocess

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



prepareDBs()
blastSearch()










#f = open('clusters')
#clusters = f.readlines()
#for cluster in clusters:
	# define dicts and lists to store cluster specific values to, and open file to write sequences of cluster to a FASTA file
#	copyDict = {}
#	seq_count = []
#	seq_list = []
#	file = open('cluster.fasta', 'w')
	
	# Count taxon specific copy numbers and store counts in a list
#	for x in set_identifiers_list:
#		count = cluster.count(x)
#		copyDict[x] = count
#		seq_count.append(count)
		#print(seq_count)
#	if sum(seq_count) > 1:
		# Store all sequence names in a list to be referred to later
#		seq = cluster.split('\t')
#		for i in seq:
#			seq_list.append(i)
		# write cluster sequences to FASTA file
		# used method below as need to differentiate between prot and nucl dicts of outgroup species. When reading the isotig (without _N or _P in name), add on the _N (only _N as dont want to read in protein seqs) by selecting name from letter 3 onwards (though need to make sure that identifiers minus the _N or _P identifiers are always 3 letters long!)
#		for i in identifiers_list:
#			if '_N' in i:
#				for x in seq:
#					if x.startswith(i.rstrip('_N')):
#						name = (i + (('').join(list(x)[3:])))
#						file.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[name] + '\n')
#	file.close()
#	if any(item.startswith('AT') for item in seq_list):
#		clusterDict = fastaDict('cluster.fasta')
#		orfs = callORF(clusterDict)

	
		




