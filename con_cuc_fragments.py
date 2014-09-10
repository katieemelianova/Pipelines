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


# For each sequence in the cucumber dict, write it to a file in FASTA format
for cucumber in X:
	blast = open('con_cuc_blastfile', 'w')
	blast.write('>' + cucumber + '\n' + X[cucumber])
	# BLASTN the cucumber sequence against the B. conchifolia transcriptome at E value cutoff of 1e-30
	subprocess.call(["blastn -query con_cuc_blastfile -db /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/CON_2.5_Isotigs.fasta -evalue 1e-30 -out con_cuc_blastoutfile -outfmt '6 qseqid sseqid qstart qend sstart send length pident evalue'"], shell=True)
	
	# open the output file of the BLASTN and check how many lines it has (how many BLAST hits)
	with open('con_cuc_blastoutfile') as myfile:
		count = sum(1 for line in myfile if line.rstrip('\n'))
		myfile.close()
	# if the cucumber sequence has hit more than one B. conchifolia sequence, continue with the script
	if count > 1:
		f = open('con_cuc_blastoutfile')
		blastout = f.readlines()
		# make an empty list to put all coordinates from all B. conchifolia hits
		total_coords = []

		for i in blastout:
			qstart = int(i.split()[2])	
			qend = int(i.split()[3])
			# turn the start and end coordinates of the hit on the cucumbers sequence into a list of numbers spanning the start and end of the hit e.g. qstart = 2 qend = 8 hit_coords = [2,3,4,5,6,7,8]
			hit_coords = list(range(qstart,qend))
			# append all the coordinates from all hits to one list
			for i in hit_coords:
				total_coords.append(i)


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

			unique = ('yes')

			# for each of the individual hit coordinates, check in the big list whether they occur only once. If not, this is assumed that they overlap with another hit. This way they are labelled eithe unique = yes or unique = no

			hit_coords = list(range(qstart,qend))
			for i in hit_coords:
				if total_coords.count(i) == 1:
					continue
				else:
					unique = ('no')
			# all unique/nonoverlapping hits are appended to a list and written to an output file
			if unique == 'yes':
				unique_list.append(con)
		
		if unique_list:
			out = open('possible_fragments', 'a')
			print(cucumber)
			out.write(cucumber + '\t')
			for i in unique_list:
				out.write(i + '\t')
			out.write('\n')

			subprocess.call(["cat possible_fragments"], shell=True)














