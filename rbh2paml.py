# rewriting original script to get kn/ks kn ks values from reciprocal best blast hits

#


import subprocess

def blast(x, y):

	subprocess.call(["blastn -query " + x + " -db " + y + " -max_target_seqs 1 -out " + namex + "_" + namey +  " -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)

	subprocess.call(["blastn -query " + y + " -db " + x + " -max_target_seqs 1 -out " + namey + "_" + namex + " -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)


def rbh(filex, filey):
	f = open(filex)
	xy = f.readlines()
	
	fi = open(filey)
	yx = fi.readlines()
	
	file = open(namex + '_' + namey + '_' + 'RBH', 'w')
	
	for i in xy:
		if i in yx:
			file.write(i)
			
			
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


#PATH = input("give me the location of the first fasta for rbh \n")
#PATH_name = input("What is the species name? \n")
x = ('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/CON_unique_fasta')
namex = ('con')

#PATH1 = input("give me the location of the second fasta for rbh \n")
#PATH_name1 = input("What is the species name? \n")
y = ('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/PLE_unique_fasta')
namey = ('ple')



blast(x, y)


f=open(namex + "_" + namey) 
xy=f.readlines()

f=open(namey + "_" + namex) 
yx=f.readlines()


# define each column, and switch first two isotig columns in one of the files, so they represent the same species as the other blast output

file_xy=open(namex + "_" + namey + '_matches', 'w')

for i in xy:
	cvp_c_scaff=i.split('\t')[0]
	cvp_p_scaff=i.split('\t')[1]
	file_xy.write(cvp_c_scaff + '\t' + cvp_p_scaff + '\n')


file_yx=open(namey + "_" + namex + '_matches', 'w')
	
for i in yx:
	pvc_p_scaff=i.split('\t')[0]
	pvc_c_scaff=i.split('\t')[1]
	file_yx.write(pvc_c_scaff + '\t' + pvc_p_scaff + '\n')

file_xy.close()
file_yx.close()

rbh(namex + "_" + namey + '_matches', namey + "_" + namex + '_matches')

subprocess.call(["rm *matches " + namex + '_' + namey + ' ' + namey + '_' + namex], shell=True)


#==========================================================================================


file = open(x)

X = {}



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

	X[C_name] = C_seq
	
file.close()


file = open(y)

Y = {}

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

	Y[C_name] = C_seq
	
file.close()


rbh2paml_outfile = open('con_ple_RBH_output', 'w')


file = open(namex + '_' + namey + '_' + 'RBH')
hits = file.readlines()

for i in hits:
	x = i.split()[0]
	y = i.split()[1]
	y = y.rstrip('\n')
	
	file = open('rbh_fasta', 'w')
	file.write('>' + x + '\n' + X[x] + '\n')
	file.write('>' + y + '\n' + Y[y] + '\n')
	
	file.close()

	file = open('rbh_fasta')
	outfile = open('rbh_fasta_longestORF', 'w')
	

	 
	
	
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
			
			
			print('>' + i + '\n' + max(FandR, key = len))
			outfile.write('>' + i + '\n' + max(FandR, key = len) + '\n')
			
			
			
	outfile.close()
	
		
			
			
	subprocess.call(["python /Volumes/BACKUP/Bioinformatics/Scripts/python_copy/BMGE-1.1/dna2pep-1.1/dna2pep.py rbh_fasta_longestORF --fasta rbh_fasta_longestORF.pep"], shell=True)
		
	subprocess.call(["linsi --quiet rbh_fasta_longestORF.pep > rbh_fasta_longestORF.pep.align"], shell=True)
		
	subprocess.call(["python /Volumes/BACKUP/Bioinformatics/Scripts/python_copy/BMGE-1.1/RevTrans-1.4/revtrans.py rbh_fasta_longestORF rbh_fasta_longestORF.pep.align rbh_fasta_longestORF.revtrans"], shell=True)
	
	
	

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
		
		rbh2paml_outfile.write(x + '\t' + y + '\t' + str(knks) + '\t' + str(kn) + '\t' + str(ks) + '\n')
	except IndexError:
		continue
		
	
	#subprocess.call(["rm rbh_fasta, rbh_fasta_longestORF, rbh_fasta_longestORF.pep, rbh_fasta_longestORF.pep.align, rbh_fasta_longestORF.revtrans"], shell=True)
	
	
	
	# put line here to remove all intermediate files so if something isnt working then it doesnt loop over the same thing
		

		# a change
	

	
	
	
	
	
	
	
