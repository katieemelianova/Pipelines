import subprocess


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



for cucumber in X:
	blast = open('con_cuc_blastfile', 'w')
	blast.write('>' + cucumber + '\n' + X[cucumber])
	subprocess.call(["blastn -query con_cuc_blastfile -db /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/CON_2.5_Isotigs.fasta -evalue 1e-30 -out con_cuc_blastoutfile -outfmt '6 qseqid sseqid qstart qend sstart send length pident evalue'"], shell=True)
	
	
	with open('con_cuc_blastoutfile') as myfile:
		count = sum(1 for line in myfile if line.rstrip('\n'))
		myfile.close()

	if count > 1:
		f = open('con_cuc_blastoutfile')
		blastout = f.readlines()

		total_coords = []

		for i in blastout:
			qstart = int(i.split()[2])	
			qend = int(i.split()[3])
			hit_coords = list(range(qstart,qend))
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

			hit_coords = list(range(qstart,qend))
			for i in hit_coords:
				if total_coords.count(i) == 1:
					continue
				else:
					unique = ('no')

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














