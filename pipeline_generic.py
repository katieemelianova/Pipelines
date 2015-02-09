import subprocess
import re
import collections



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

def mcl(file, inflation):
	# Calls MCL function to create index files necessary for MCL clustering from the concatenated FASTA file of all species
	subprocess.call(["mcxload -abc " + file + " --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o " + file.rstrip('.abc') + ".mci -write-tab " + file.rstrip('.abc') + ".tab"], shell=True)
	# Calls MCL. The value following '-I ' is the inflation parameter, larger -I yields fewer clusters. 
	subprocess.call(["mcl " + file.rstrip('.abc') + ".mci -I " + str(inflation)], shell=True)

	# Uses index files to create the output file using sequence headers from all species. If you changed the inflation parameter in the previous stage, change the two occurences of it in this command also. e.g. -I 1.5 == .mci.I15
	subprocess.call(["mcxdump -icl out." + file.rstrip('.abc') + ".mci.I30 -tabr "+ file.rstrip('.abc') + ".tab -o dump." + file.rstrip('.abc') + ".mci.I30"], shell=True)

def blast(file):
	# makes a BLAST directory and uses discontinuous megablast to BLAST concatenated file of all species against itself [NEEDS TO BE GENERALIZED]
	subprocess.call(["makeblastdb -in /home/kemelianova/ncbi-blast-2.2.27+/db/alltranscriptomes2blast -dbtype nucl"], shell=True)
	subprocess.call(["blastn -task dc-megablast -query /home/kemelianova/ncbi-blast-2.2.27+/db/alltranscriptomes2blast -db /home/kemelianova/ncbi-blast-2.2.27+/db/alltranscriptomes2blast -out alltranscriptomes2mcl.blastout -evalue 1e-70 -outfmt '6 qseqid sseqid evalue length pident'"], shell=True)
	subprocess.call(["mv alltranscriptomes2mcl.blastout alltranscriptomes2mcl.abc"], shell=True)
	subprocess.call(["rm /home/kemelianova/ncbi-blast-2.2.27+/db/alltranscriptomes2blast*"], shell=True)

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

def codons(seq):
	# Finds longest ORF (here defined as longest sequence with no stop codons in any frame)
	try:
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
	except ValueError:
		print('I dont think there are any sequences here')	


def dists():
	# Calls mothur, estimates pairwise distances, and then records the pairwise distances for each sequence to every other sequence.
	subprocess.call(['mothur "#dist.seqs(fasta=cluster.orfs.transalign);"'], shell=True)
	f = open('cluster.orfs.dist')
	dists = f.readlines()

	seqs = []

	for i in dists:
		first_seq = i.split()[0]
		second_seq = i.split()[1]
		seqs.append(first_seq)
		seqs.append(second_seq)

	seqs = set(seqs)

	seqDists = {}

	for i in seqs:
		currentSeqDists = []
		for line in dists:
			if i in line:
				currentSeqDist = float(line.split()[2])
				currentSeqDists.append(currentSeqDist)
		meanDist = sum(currentSeqDists)/len(currentSeqDists)
		seqDists[i] = meanDist
	return(seqDists)

def consensus():
	# Calls consensus.R, parses the output file to get the percentage of bases which have consensus under the threshold set in script consensus.R
	subprocess.call(["Rscript consensus.R"], shell=True)			
	f = open('alignment_consensus')
	consensus = f.readlines()
	consensus_bases = 0
	total_bases = 0
	for i in consensus:
		i = i.rstrip('\n')
		total_bases+=1
		if i != 'NA':
			if i != '-':
				consensus_bases+=1
	percent_consensus_bases = (float(consensus_bases*100)/float(total_bases))
	return(percent_consensus_bases)


# different options for alignment of protein sequences and translational aligner
def dna2pep():
	subprocess.call(["python dna2pep-1.1/dna2pep.py --outformat fasta cluster_longest_orfs > cluster.orfs.pep"], shell=True)
def mafft():
	subprocess.call(["mafft --quiet cluster.orfs.pep > cluster.orfs.pep.aln"], shell=True)
def tcoffee():
	subprocess.call(["t_coffee cluster.orfs.pep -output=fasta_aln"], shell=True)
	subprocess.call(["mv cluster.orfs.fasta_aln cluster.orfs.pep.aln"], shell=True)
def revtrans():
	subprocess.call(["python RevTrans-1.4/revtrans.py cluster_longest_orfs cluster.orfs.pep.aln cluster.orfs.transalign"], shell=True)
def pal2nal():
	subprocess.call(["perl pal2nal.v14/pal2nal.pl -output fasta cluster.orfs.pep.aln cluster_longest_orfs > cluster.orfs.transalign"], shell=True)



def del_intermediates():
	subprocess.call(['rm alignment_consensus'], shell=True)
	subprocess.call(['rm *cluster.orf.raxml*'], shell=True)
	subprocess.call(['rm paml_out'], shell=True)
	subprocess.call(['rm mothur.*'], shell=True)
	subprocess.call(['rm erraafile.*'], shell=True)
	subprocess.call(['rm errnucfile.*'], shell=True)
	subprocess.call(['rm RAxML_bestTree.cluster.orf.raxml'], shell=True)
	subprocess.call(['rm RAxML_log.cluster.orf.raxml'], shell=True)
	subprocess.call(['rm RAxML_result.cluster.orf.raxml'], shell=True)
	subprocess.call(['rm RAxML_info.cluster.orf.raxml'], shell=True)
	subprocess.call(['rm RAxML_parsimonyTree.cluster.orf.raxml'], shell=True)
	              
      

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
	path=i.split()[1]	
	filename = open(path)
	dictionary_seq_dictionaries['{0}'.format(identifier)]=dictionary_seq(identifier, filename)
	filename.close()

# write contents of dict of dicts to a multifasta file, filtering by sequence length
all_trans = open('/home/kemelianova/ncbi-blast-2.2.27+/db/alltranscriptomes2blast', 'w')
for i in identifiers_list:
	for x in (dictionary_seq_dictionaries[i]):
		if int(len((dictionary_seq_dictionaries[i])[x])) > 200:
			all_trans.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')



blast('/home/kemelianova/ncbi-blast-2.2.27+/db/alltranscriptomes2blast')
mcl('alltranscriptomes2mcl.abc', 3)


# Open the outfile to write results to
outfile = open('pipline_generic_outfile', 'w')

f = open('dump.alltranscriptomes2mcl.mci.I30')
dump = f.readlines()
for cluster in dump:
	# define dicts and lists to store cluster specific values to, and open file to write sequences of cluster to a FASTA file
	copyDict = {}
	seq_count = []
	seq_list = []
	file = open('cluster_seqs', 'w')
	
	# Count taxon specific copy numbers and store counts in a list
	for x in identifiers_list:
		count = cluster.count(x)
		copyDict[x] = count
		seq_count.append(count)
	# Store all sequence names in a list to be referred to later
	seq = cluster.split()
	for i in seq:
		seq_list.append(i)
	# write cluster sequences to FASTA file
	for i in identifiers_list:
		for x in seq:
			if x.startswith(i):
				file.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')

	file.close()

	# Limit further analyses to clusters with more than one sequence
	if sum(seq_count) > 1:

		file = open('cluster_seqs')
		out = open('cluster_longest_orfs', 'w')

		# for each cluster, get the longest ORF and write to a new file

		L = {}

		for line in file:
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


		for i in L:
			s = L[i]
			rcs = rc(s)
			stop = ('TAG','TAA','TGA')
			# check whether the sequence has no stop codons anywhere
			if any(i not in s for i in stop) or any(i not in rcs for i in stop):
				out.write('>' + i + '\n' + s + '\n')
			else:
				FandR = []
				F = codons(s)
				R = codons(rc(s))
				FandR.append(F)
				FandR.append(R)
				out.write('>' + i + '\n' + max(FandR, key = len) + '\n')
				
		# translate, align and translationall align all sequences

		dna2pep()
		mafft()
		pal2nal()

		f = open('cluster.orfs.transalign')
		test = f.read()
		if test: 
			f.close()
		
			# Make a ML tree using RAxML if cluster size exceeds 3
			if sum(seq_count) > 3:

				subprocess.call(["./standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRCAT -p 12345 -s cluster.orfs.transalign -n cluster.orf.raxml"], shell=True)
			
			# Write taxon copy numbers to oufile
			for i in identifiers_list:
				outfile.write(str(copyDict[i]) + '\t')

			# Write % consensus bases to outfile
			consensus_bases = consensus()
			outfile.write(str(consensus_bases) + '\t')

			# Write all pairwise distances to outfile
			distDict = dists()
			outfile.write(('| start dists | ') + '\t')
			for i in distDict:
				dist = str(distDict[i])
				outfile.write(dist + '\t')
			outfile.write(('| end dists | ') + '\t')				
			

			# Call codeml
			subprocess.call(['codeml'], shell=True)

			f = open('paml_out')
			paml = f.read()

			# Use try statement to prevent stalling due to no codeml result
			try:
				# Use split to isolate matrix output by codeml from other output and remove brackets and values within them (these contain kN and kS values)
				knks_list = []
				output = paml.split('comparison.)\n\n')
				output = (output[1])
				output = output.split('TREE')
				output = output[0]
				knks = re.sub(r'\([^)]*\)', '', output)
				
				# Add all values from the matrix to a list, excluding sequence names and values returned with a -1.0000 value
				for i in knks.split():
					if not any(x in i for x in identifiers_list):
						if i != '-1.0000':
							knks_list.append(i)
				
				# Reformat bracket placement for ease of splitting
				output = output.replace(')-', ') -')
				kn_list = []
				ks_list = []
				# Write kN and kS values to separate lists
				for i in output.split():
					if '(' in i:
						kn = i.lstrip('(')
						if kn != '-1.0000':
							kn_list.append(kn)
					elif ')' in i: 
						ks = i.rstrip(')')
						if ks != '-1.0000':
							ks_list.append(ks)
				# Write kN/kS, kN and kS to outfile
				outfile.write('|knks start|' + '\t')
				for i in knks_list:
					outfile.write(str(i) + '\t')
				outfile.write('|knks end|' + '\t')
				outfile.write('|kn start|' + '\t')
				for i in kn_list:
					outfile.write(str(i) + '\t')
				outfile.write('|kn end|' + '\t')
				outfile.write('|ks start|' + '\t')
				for i in ks_list:
					outfile.write(str(i) + '\t')
				outfile.write('|ks end|' + '\t')
				outfile.write('|sequence name start|' + '\t')
				for i in seq_list:
					outfile.write(i + '\t')
				outfile.write('|sequence name end|' + '\t')
				outfile.write('\n')
			except IndexError:
					print('no paml result')
					outfile.write('\n')

			# Delete all intermediate files
			del_intermediates()
		
