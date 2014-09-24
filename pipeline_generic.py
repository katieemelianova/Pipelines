import subprocess
import re
import collections



def dictionary_seq(identifier, filename):
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

def mcl(file):
	subprocess.call(["mcxload -abc " + file + " --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o " + file.rstrip('.abc') + ".mci -write-tab " + file.rstrip('.abc') + ".tab"], shell=True)
	
	subprocess.call(["mcl " + file.rstrip('.abc') + ".mci -I 3"], shell=True)

	subprocess.call(["mcxdump -icl out." + file.rstrip('.abc') + ".mci.I30 -tabr "+ file.rstrip('.abc') + ".tab -o dump." + file.rstrip('.abc') + ".mci.I30"], shell=True)

def blast(file):
	subprocess.call(["makeblastdb -in /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/alltranscriptomes2blast -dbtype nucl"], shell=True)
	subprocess.call(["blastn -query /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/alltranscriptomes2blast -db /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/alltranscriptomes2blast -out alltranscriptomes2mcl.abc -outfmt '6 qseqid sseqid evalue'"], shell=True)
	subprocess.call(["rm /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/alltranscriptomes2blast*"], shell=True)

def rc(seq): 
	
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'} 
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


dictionary_seq_dictionaries={}
identifiers_list=[]

f = open('control_testfile')
testfile = f.readlines()


for i in testfile:
	identifier=i.split()[0]
	identifier = identifier.upper()
	identifiers_list.append(identifier)
	path=i.split()[1]	
	filename = open(path)
	dictionary_seq_dictionaries['{0}'.format(identifier)]=dictionary_seq(identifier, filename)
	filename.close()

all_trans = open('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/alltranscriptomes2blast', 'w')
for i in identifiers_list:
	for x in (dictionary_seq_dictionaries[i]):
		all_trans.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')



blast('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/alltranscriptomes2blast')
mcl('alltranscriptomes2mcl.abc')

outfile = open('pipline_generic_outfile', 'w')

f = open('dump.alltranscriptomes2mcl.mci.I30')
dump = f.readlines()
for cluster in dump:
	file = open('cluster_seqs', 'w')
	for x in identifiers_list:
		count = i.count(x)
		outfile.write(str(count) + '\t')
	seq = cluster.split()
	for i in identifiers_list:
		for x in seq:
			if x.startswith(i):
				file.write('>' + x + '\n' + (dictionary_seq_dictionaries[i])[x] + '\n')

	file.close()



	file = open('cluster_seqs')
	out = open('cluster_longest_orfs', 'w')

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
			


	subprocess.call(["python dna2pep-1.1/dna2pep.py cluster_longest_orfs --fasta cluster.orfs.pep"], shell=True)
		
	subprocess.call(["linsi --quiet cluster.orfs.pep > cluster.orfs.pep.align"], shell=True)
		
	subprocess.call(["python RevTrans-1.4/revtrans.py cluster_longest_orfs cluster.orfs.pep.align cluster.orfs.revtrans"], shell=True)

	subprocess.call(["/Volumes/BACKUP/Bioinformatics/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -T 2 -m GTRCAT -p 12345 -s cluster.orfs.revtrans -n cluster.orf.raxml"], shell=True)

	subprocess.call(["Rscript consensus.R"], shell=True)
		
	f = open('alignment_consensus')
	consensus = f.readlines()
	
	consensus_bases = 0
	total_bases = 0
	
	for i in consensus:
		total_bases+=1
		if i != 'NA':
			if i != '-':
				consensus_bases+=1
	percent_consensus_bases = (float(consensus_bases*100)/float(total_bases))

	subprocess.call(["cat cluster.orfs.revtrans"], shell=True)


	subprocess.call(['codeml'], shell=True)

	f = open('paml_out')


	paml = f.read()

	try:
		knks_list = []
		kn_ks_list = []
		output = paml.split('comparison.)\n\n')
		output = (output[1])
		output = output.split('TREE')
		output = output[0]
		knks = re.sub(r'\([^)]*\)', '', output)
		for i in knks.split():
			if not any(x in i for x in identifiers_list):
				if i != '-1.0000':
					knks_list.append(i)
		
		output = output.replace(')-', ') -')
		kn_ks_list = []
		for i in output.split():
			if ')' in i or '(' in i:
				kn_ks = i.lstrip('(')
				kn_ks = kn_ks.rstrip(')')
				if kn_ks != '-1.0000':
					kn_ks_list.append(kn_ks)

	except IndexError:
			print('no paml result')


	subprocess.call(['rm *cluster.orf.raxml*'], shell=True)

