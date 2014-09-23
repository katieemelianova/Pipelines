#!/usr/bin/python






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


#path = input(' what is the full path name of the fasta file?')
file = ('longest_ORF_KE_problemtest.txt')

f = open(file)


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


for i in L:
	s = L[i]
	rcs = rc(s)
	stop = ('TAG','TAA','TGA')
	# if no stop codons in the sequence, print the sequence
	if any(i not in s for i in stop) or any(i not in rcs for i in stop):
		print('>' + i + '\n' + s + '\n')
	else:
		FandR = []
		F = codons(s)
		R = codons(rc(s))
		FandR.append(F)
		FandR.append(R)
		print('>' + i + '\n' + max(FandR, key = len))
		

	
	
	
	
	
	
		
	
	
	
	
	
	
	
	
