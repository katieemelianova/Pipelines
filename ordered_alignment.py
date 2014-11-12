#ordered_alignment.py

import subprocess
import collections


PATH = input(' give me full path of the Multi-FASTA file you want to align ')
name = (PATH)

subprocess.call(["sizeseq -sequences "+ name + " -outseq " + name + "_sizeseq -descending Y"], shell=True)

file = open(name + "_sizeseq")

C_name2seq = collections.OrderedDict()

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
	
	C_name2seq[C_name] = C_seq



count = 0
for i in C_name2seq:
	count+=1
	if count <= 2:
		file = open('seq2align', 'a')
		file.write('>' + i + '\n' + C_name2seq[i] + '\n')
		file.close()
		if count == 2:
			subprocess.call(["mafft seq2align > seq.mafft"], shell=True)
	elif count > 2:
		file = open('seq2align', 'w')
		file.write('>' + i + '\n' + C_name2seq[i] + '\n')
		file.close()
		subprocess.call(["mafft --add seq2align --reorder seq.mafft > seq.mafft_temp"], shell=True)
		subprocess.call(["mv seq.mafft_temp seq.mafft"], shell=True)





















