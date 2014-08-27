# modification of the rbh2paml.py script. This version has the same method of alignment and running through PAML, however it finds paralogs instead of orthologs.
#Method of paralog identification is: blast transcriptome against itself at 1e-40, with maxtargetseqs = 2. Filter out all hits which are a sequence against itself and also reciprocal hits.
# Filter out all hits which have a hit length of less than 200bp.
# Continue with script as in rbh2paml.py

import subprocess

def blast(x):

	subprocess.call(["blastn -query " + x + " -db " + x + " -max_target_seqs 2 -evalue 1e-40 -out " + namex +  "_paralog_blast -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)

def parse():
	f = open(namex +  "_paralog_blast")
	blast = f.readlines()
	outfile = open(namex + '_paralogs', 'w')
	
	hit_list = []
	
	for i in blast:
		q = i.split()[0]
		d = i.split()[1]
		qd = (q + '/t' + d)
		hit_list.append(qd)
		dq = (d + '/t' + q)
		if dq not in hit_list:
			l = i.split()[2]
			pi = i.split()[3]
			if q!=d:
				if int(l) > 200:
					outfile.write(i)
			
def rc(seq): 
	
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-', 'W':'W', 'Y':'Y', 'M':'M', 'S':'S', 'K':'K'} 
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
		

		
		

		

		
		
	
		
		
		
		
	



#===================================================


PATH = input("give me the location of the fasta for paralog identification \n")
PATH_name = input("What is the species name? \n")
x = (PATH)
namex = (PATH_name)




blast(x)
parse()




#==========================================================================================


file = open(x)

dict = {}



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

	dict[C_name] = C_seq
	
file.close()



paralog2paml_outfile = open(namex + '_' +'paralog_output', 'w')


file = open(namex + '_paralogs')
hits = file.readlines()

for i in hits:
	x = i.split()[0]
	y = i.split()[1]
	y = y.rstrip('\n')

	
	file = open('paralog_fasta', 'w')
	file.write('>' + x + '\n' + dict[x] + '\n')
	file.write('>' + y + '\n' + dict[y] + '\n')
	
	file.close()

	file = open('paralog_fasta')
	outfile = open('paralog_fasta_longestORF', 'w')
	

	 
	
	
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
	
	file.close()

	
	for i in L:
		s = L[i]
		rcs = rc(s)
		stop = ('TAG','TAA','TGA')
		if any(i not in s for i in stop) or any(i not in rcs for i in stop):
			continue
			#print('>' + i + '\n' + s + '\n')
		else:
			FandR = []
			F = codons(s)
			R = codons(rc(s))
			FandR.append(F)
			FandR.append(R)
			
			
			
			outfile.write('>' + i + '\n' + max(FandR, key = len) + '\n')
			
			
			
	outfile.close()
	
		
			
			
	subprocess.call(["python dna2pep-1.1/dna2pep.py paralog_fasta_longestORF --fasta paralog_fasta_longestORF.pep"], shell=True)
		
	subprocess.call(["linsi --quiet paralog_fasta_longestORF.pep > paralog_fasta_longestORF.pep.align"], shell=True)
		
	subprocess.call(["python RevTrans-1.4/revtrans.py paralog_fasta_longestORF paralog_fasta_longestORF.pep.align paralog_fasta_longestORF.revtrans"], shell=True)
	
	
	

	subprocess.call(['codeml'], shell=True)
	
	try:
		f = open ('paml_out')
		lines = f.readlines()
		line_list = []
		for i in lines:
			line_list.append(i)
		last = line_list[-1]
		knks = last.split()[1]
		kn = last.split()[2]
		kn = kn.lstrip('(')
		ks = last.split()[3]
		ks = ks.rstrip(')')
		
		print(knks)
		print(kn)
		print(ks)
		
		paralog2paml_outfile.write(x + '\t' + y + '\t' + str(knks) + '\t' + str(kn) + '\t' + str(ks) + '\n')
	except IndexError:
		continue
		
	
	subprocess.call(["rm paralog_fasta, paralog_fasta_longestORF, paralog_fasta_longestORF.pep, paralog_fasta_longestORF.pep.align, paralog_fasta_longestORF.revtrans"], shell=True)
	

	
	
	
	
	
	
	
