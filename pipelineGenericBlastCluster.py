import subprocess

def blastn(file):
	# makes a BLAST directory and uses discontinuous megablast to BLAST concatenated file of all species against itself [NEEDS TO BE GENERALIZED]
	subprocess.call(["makeblastdb -in alltranscriptomes2blastn -dbtype nucl"], shell=True)
	subprocess.call(["blastn -task dc-megablast -query alltranscriptomes2blastn -db alltranscriptomes2blastn -out alltranscriptomes.blastnout -evalue 1e-40 -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)

def blastx(nuc, prot):
	# makes a BLAST directory and uses discontinuous megablast to BLAST concatenated file of all species against itself [NEEDS TO BE GENERALIZED]
	subprocess.call(["makeblastdb -in alltranscriptomes2blastx -dbtype prot"], shell=True)
	subprocess.call(["blastx -query alltranscriptomes2blastn -db alltranscriptomes2blastx -out alltranscriptomes.blastxout -evalue 1e-40 -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)
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


# create a dict to hold all the dicts for all species
dictionary_seq_dictionaries={}
# Create an empty list to hold all species identifiers to refer to later
identifiers_list=[]

f = open('pipeline_control')
testfile = f.readlines()
# Read in contents of control file, creating a dict for each species, and adding these dicts to another dict
for i in testfile:
	identifier=i.split()[0]
	identifier = identifier.upper()
	identifiers_list.append(identifier)
	identifier = identifier.split('_')[0]
	path=i.split()[1]	
	filename = open(path)
	dictionary_seq_dictionaries['{0}'.format(identifier)]=dictionary_seq(identifier, filename)
	filename.close()

# write contents of dict of dicts to a multifasta file, filtering by sequence length
all_transN = open('alltranscriptomes2blastn', 'w')
for i in identifiers_list:
	if '_N' in i:
		i = i.split('_')[0]
		for x in (dictionary_seq_dictionaries[i]):
			if int(len((dictionary_seq_dictionaries[i])[x])) > 100:
				all_transN.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')

# write contents of dict of dicts to a multifasta file, filtering by sequence length
all_transP = open('alltranscriptomes2blastx', 'w')
for i in identifiers_list:
	if '_OP' in i:
		i = i.split('_')[0]
		for x in (dictionary_seq_dictionaries[i]):
			if int(len((dictionary_seq_dictionaries[i])[x])) > 100:
				all_transP.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')

blastn('alltranscriptomes2blastn')
blastx('alltranscriptomes2blastn', 'alltranscriptomes2blastx')


#clustering step

#subprocess.call(["cat alltranscriptomes.blastxout alltranscriptomes.blastnout > alltranscriptomes.blastout"], shell=True)

#hitLists = []
#f = open('alltranscriptomes.blastout')
#blast = f.readlines()
#for i in blast:
#	q = i.split()[0]
#	d = i.split()[1]
#	l = i.split()[2]
#	p = i.split()[3]
#	e = i.split()[4]
	# limit to searches that arent self-blasted, have a length of > 150 and have > 40% identity. Add all match pairs with this description to a list of lists
#	if q != d:
#		if int(l) > 150:
#			if float(p) > 40:
#				hitPair = [q,d]
#				hitLists.append(hitPair)
# go through each pair and either merge with a another group, if one or more sequences are shared, or add to the list of groups. Write finished groups to a list
#merged = []
#outFile = open('clusters', 'w')
#for lists in [set(d) for d in hitLists]:
#	for m in merged:
#		if not m.isdisjoint(lists):
#			m.update(lists)
#			break
#	else:
#		merged.append(lists)
#for i in merged:
#	for s in i:
#		outFile.write(s.rstrip('\n') + '\t')
#	outFile.write('\n')

#f = open('clusters')
#clusters = f.readlines()
#for cluster in clusters:
	# define dicts and lists to store cluster specific values to, and open file to write sequences of cluster to a FASTA file
#	copyDict = {}
#	seq_count = []
#	seq_list = []
#	file = open('cluster.fasta', 'w')
	
	# Count taxon specific copy numbers and store counts in a list
#	for x in identifiers_list:
#		count = cluster.count(x)
#		copyDict[x] = count
#		seq_count.append(count)
#	if sum(seq_count) > 1:
		# Store all sequence names in a list to be referred to later
#		seq = cluster.split('\t')
#		for i in seq:
#			seq_list.append(i)
		# write cluster sequences to FASTA file
#		for i in identifiers_list:
#			for x in seq:
#				if x.startswith(i):
#					file.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')

#	file.close()

#	clusterDict = fastaDict('cluster_seqs')
#	out = open('clusterLongestOrfs', 'w')
#	for i in L:
#		s = L[i]
#		rcs = rc(s)
#		stop = ('TAG','TAA','TGA')
		# check whether the sequence has no stop codons anywhere
#		if any(i not in s for i in stop) or any(i not in rcs for i in stop):
#			out.write('>' + i + '\n' + s + '\n')
#		else:
#			FandR = []
#			F = codons(s)
#			R = codons(rc(s))
#			FandR.append(F)
#			FandR.append(R)
#			out.write('>' + i + '\n' + max(FandR, key = len) + '\n')

#	outfile = open('clusterLongestOrfs', 'w')
#	orfDict = fastaDict('clusterLongestOrfs')
#	for i in orfDict:
#		outfile.write('>' + i + '\n' + orfDict[i] + '\n')
	




