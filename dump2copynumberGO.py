import subprocess


file = open('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/CON_unique_fasta')
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

file = open('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/PLE_unique_fasta')
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

file = open('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/VEN_unique_fasta')
Z = {}
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

	Z[C_name] = C_seq
	
file.close()

f = open('/Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/TAIR10_pep_annotations')
AT = {}
ath = f.readlines()
for i in ath:
	if i.startswith('>'):
		at = (i.split('|')[0]).lstrip('>')
		at = at.rstrip(' ')
		go = i.split('|')[2]
		
		AT[at] = go



f = open('dump.conpleven_selfie.mci.I30')
dump = f.readlines()


outfile = open('conpleven_copynumberGO', 'w')

for i in dump:
	fasta2blast = open('dump2cpGO_fasta', 'w')
	out = open('dump_copynumerGO', 'w')
	con_count = i.count('CON_')
	ple_count = i.count('PLE_')
	ven_count = i.count('VEN_')
	members = i.split()
	for member in members:
		if member.startswith('CON_'):
			member  = member.lstrip('CON_')	
			fasta2blast.write('>CON' + member + '\n' + X[member.lstrip('CON_')] + '\n')
			
		elif member.startswith('PLE_'):
			member = member.lstrip('PLE_')
			fasta2blast.write('>PLE' + member + '\n' + Y[member.lstrip('PLE_')] + '\n')

		elif member.startswith('VEN_'):
			member = member.lstrip('VEN_')
			fasta2blast.write('>VEN' + member + '\n' + Z[member.lstrip('VEN_')] + '\n')
	fasta2blast.close()
	#fasta2blast = open('dump2cpGO_fasta')
	#fasta = fasta2blast.read()
	#print(fasta)
	#print('\n')
	subprocess.call(["blastx -query dump2cpGO_fasta -db /Volumes/BACKUP/Bioinformatics/ncbi-blast-2.2.27+/db/TAIR10_pep_annotations -evalue 1e-60 -out dump_blast -outfmt '6 qseqid sseqid pident evalue'"], shell=True)
		
		
	f = open('dump_blast')
	blast = f.readlines()
	at_GO_list = []
	for i in blast:
		at = i.split()[1]
		at_GO_list.append(AT[at])
	
	outfile.write(str(con_count) + '\t' + str(ple_count) + '\t' + str(ven_count) + ' | ')
	for i in (set(at_GO_list)):
		outfile.write(i + '\t')
	outfile.write('|' + '\t')
	for member in members:
		outfile.write(member + '\t') 
	outfile.write('\n')	
		
			
		
			
			
		
	


	
	