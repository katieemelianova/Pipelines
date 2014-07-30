#!/usr/bin/python

# This script uses two FASTA files to conduct a reciprocal best BLAST hit analysis to identify putative orthologs. The longest ORF is extracted from each ortholog and used to create a translational alignment. This alignment is used as input for PAML, giving an output of ortholog pair IDs, kN/kS, kN, and kS values.

# the packages required to run this script are: EMBOSS, BLAST, dna2pep.py, revtrans.py, mafft and PAML. The PAML (package codeml) control file uses a one ratio model (model = 0).

# Before starting, move the two FASTA files you are using in the analysis into your BLAST db folder, and use the makeblastdb command to make blast databases for the files.

# EMBOSS, BLAST and mafft can be called from any directory. Ensure dna2pep and Revtrans are in the same directory (unless you reference the full path name in the script). PAML can be called from the working directory as long as the control file (codeml.ctl) is located in the same directory.


import re
import collections
import subprocess

def blast(x, y):
	# make two reciprocal BLAST searches, with a maximum target sequence of 1 and no e-value threshold.
	subprocess.call(["blastn -query " + x + " -db " + y + " -max_target_seqs 1 -out " + namex + "_" + namey +  " -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)

	subprocess.call(["blastn -query " + y + " -db " + x + " -max_target_seqs 1 -out " + namey + "_" + namex + " -outfmt '6 qseqid sseqid length pident evalue'"], shell=True)

def rbh(filex, filey):
	
	# Use the two files with the same species/column identities to find where gene pairs have reciprocal best BLAST hits.
	
	f = open(filex)
	xy = f.readlines()
	
	fi = open(filey)
	yx = fi.readlines()
	
	file = open(namex + '_' + namey + '_' + 'RBH', 'w')
	
	for i in xy:
		if i in yx:
			file.write(i)
			

def getorf():
	# Call getorf, looking for nucleic sequences between start and stop codons
	subprocess.call(["getorf -sequence rbh_fasta -find 3 -auto -outseq " + 'orfs'], shell=True)


def sizeseq():
	# Call sizeseq, sorting all ORFs found in getorf output by descending length
	subprocess.call(["sizeseq -sequences orfs -descending Y -outseq orfs_sorted"], shell=True)


def format():
	# Using an ordered dict to maintain descending size order of sequences, change the output of sizeseq to a format easier to find and remove duplicates.
	file = open('orfs_sorted')	
		
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
			C_seq = C_seq.rstrip('\n')
			C_seq = C_seq.upper()
	
		C_name2seq[C_name] = C_seq
			
	
		f = open('for_dup_removal', 'w')
	
	
	for k, v in C_name2seq.items():
		f.write('>' + k + '  ' + v + '\n')

def remove_duplicates():
	# Duplicates will be found due to multiple ORFs being found in each sequence. 
	# These are removed by looking at each sequence in size descending order and writing them to a list. Each sequence name checked is only taken if the name does not appear in the list. This way, the longest ORF from each sequence is retained and the rest discarded.
	f = open('for_dup_removal')
	fasta  = f.read()
	no_underscore = re.sub(r'\_.*?\ ', '', fasta)
	item = (no_underscore.split('\n'))
	lists = []
	file = open('longest_orfs', 'w')
	item = list(filter(None, item))
	for i in item:
		name = (i.split()[0])	
		seq = 	(i.split()[1])
		for i in name:
			while name not in lists:
				lists.append(name)
				file.write(name + '\n' + seq + '\n')

#=====================Reciprocal Best BLAST hits Start===============================================

# Enter the full path of the FASTA files you want to use, and give them short memorable identifiers.

PATH = input("give me the location of the first fasta for rbh \n")
PATH_name = input("What is the species name? \n")
x = (PATH)
namex = (PATH_name)
PATH1 = input("give me the location of the second fasta for rbh \n")
PATH_name1 = input("What is the species name? \n")
y = (PATH1)
namey = (PATH_name1)

blast(x, y)

f=open(namex + "_" + namey) 
xy=f.readlines()

f=open(namey + "_" + namex) 
yx=f.readlines()

# switch columns in one of the files so they represent the same species as the reciprocal BLAST output file.

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

#===========================Reciprocal Best BLAST hits End=======================================


#===========================Retrieve longest ORF Start=======================================

# Make dicts of both original FASTA files {name:sequence}

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

# For each ortholog pair, use the name to retrieve the sequence from appropriate dict and create FASTA file of ortholog pair.

for i in hits:
	x = i.split()[0]
	y = i.split()[1]
	y = y.rstrip('\n')
	file = open('rbh_fasta', 'w')
	file.write('>' + x + '\n' + X[x] + '\n')
	file.write('>' + y + '\n' + Y[y] + '\n')
	
	file.close()
	getorf()
	sizeseq()
	format()
	remove_duplicates()
	
#===========================Retrieve longest ORF End=======================================

#===========================Alignment and codeml Analysis Start=================================
		
	# Because getorf uses start and stop codons to identify ORFs, sequences with no start or stop codons are not processed, preventing downstream analyses e.g. alignment.
	# The 'longest_orfs' file is opened to check it has two sequences in it, giving an error message if this is not the case. The ortholog pairs producing such error messages are discarded, but the script does not break.
	file = open('longest_orfs')
	orfs = file.read()
	if int(orfs.count('>')) == 2:
		file.close()
			
		# ortholog pairs are translated using dna2pep.py, peptide sequences aligned using mafft, and revtrans.py is used to translationally align the nucleotide sequences using the peptide alignment as a guide.
		
		subprocess.call(["python dna2pep-1.1/dna2pep.py longest_orfs --fasta rbh_fasta_longestORF.pep"], shell=True)
			
		subprocess.call(["linsi --quiet rbh_fasta_longestORF.pep > rbh_fasta_longestORF.pep.align"], shell=True)
			
		subprocess.call(["python RevTrans-1.4/revtrans.py longest_orfs rbh_fasta_longestORF.pep.align rbh_fasta_longestORF.revtrans"], shell=True)
		
		# codeml package from PAML is called.
		
		subprocess.call(['codeml'], shell=True)
		
		# the output file of codeml is used to extract the kN/kS ratio, and kN and kS values. 
		# A try/except statement is used in case there is a problem with the codeml analysis and there is no output to stop the script from breaking.
		
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

			# Write to output file sequence IDs of ortholog pair, kN/kS ratio, and kN and kS values.
			rbh2paml_outfile.write(x + '\t' + y + '\t' + str(knks) + '\t' + str(kn) + '\t' + str(ks) + '\n')
		except IndexError:
			continue
	else:
		print('\n' + 'At least one of your gene pairs does not have a long enough ORF. This gene pair has been discarded.' + '\n')
		
	#===========================Alignment and codeml Analysis End=================================
	
	# Remove all intermediate files generated by script.	
	subprocess.call(["rm orfs orfs_sorted for_dup_removal longest_orfs rbh_fasta longest_orfs rbh_fasta_longestORF.pep rbh_fasta_longestORF.pep.align rbh_fasta_longestORF.revtrans paml_out"], shell=True)