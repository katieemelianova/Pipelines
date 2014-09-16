#!usr/bin/python

# This script (cuurently) uses the B. conchifolia unique transcriptome (alternative transcripts removed) and the cucumber (primrary transcripts only) transcriptome.

# The purpose of the script is to identify cucumber transcripts which have two ranscripts matching to it in non-overlapping regions.

import subprocess


# Create a dict of the cucumber sequences
c = open('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/Csativus_122_transcript_primaryTranscriptOnly.fa')

X = {}
for line in c:
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
c.close()

out = open('possible_fragments', 'w')


# For each sequence in the cucumber dict, write it to a file in FASTA format
for cucumber in X:
	blast = open('con_cuc_blastfile', 'w')
	blast.write('>' + cucumber + '\n' + X[cucumber])
	# BLASTN the cucumber sequence against the B. conchifolia transcriptome at E value cutoff of 1e-30
	subprocess.call(["blastn -query con_cuc_blastfile -db /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/CON_unique_fasta -evalue 1e-30 -out con_cuc_blastoutfile -outfmt '6 qseqid sseqid qstart qend sstart send length pident evalue'"], shell=True)

	f = open('con_cuc_blastoutfile')
	blastout = f.readlines()
	# make an empty list to put all coordinates from all B. conchifolia hits
	coords_lists = []
	for i in blastout:
		con = i.split()[1]
		cuc  = i.split()[0]
		qstart = int(i.split()[2])	
		qend = int(i.split()[3])
		# turn the start and end coordinates of the hit on the cucumbers sequence into a list of numbers spanning the start and end of the hit e.g. qstart = 2 qend = 8 hit_coords = [2,3,4,5,6,7,8]
		# append all the coordinates from all hits to one list
		coords_lists.append(list(range(qstart,qend)))


	unique_list = []
	for i in blastout:
		cuc = i.split()[0]
		con = i.split()[1]
		qstart = int(i.split()[2])	
		qend = int(i.split()[3])
		sstart = i.split()[4]
		send = i.split()[5]
		pident = i.split()[6]
		evalue = i.split()[7]
		
		# set hit_coords variable as list of numbers spanning length of hit in this loop. For each hit coord stored in coords_lists, get transcripts which do not overlap with this one. Repeat for each hit in the blast outfile.
		hit_coords = list(range(qstart,qend))
		for i in coords_lists:
			if not ((set(i)) & set(hit_coords)):
				unique_list.append(con)
			else:
				continue
	
	# Get all unique hits ientified and write them to a file, unless there is only one hit in the file. This can happen because one transcript hits a cucumber sequence in two non-overlapping areas,. Also removes single sequence hits.

	if len(set(unique_list)) > 1:
		out.write(cuc + '\t')
		for i in (set(unique_list)):
			out.write(i + '\t')
		out.write('\n')
