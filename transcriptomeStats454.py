import subprocess
import collections


# a script for generating some transcriptome metrics
# N50, % seqs < 200bp, % seqs < 1kb, % seqs > 1kb, min seq length, max seq length


def n50(name):
	subprocess.call(["sizeseq -sequences "+ name + " -outseq " + name + "_sizeseq -descending Y"], shell=True)

	files = open(name + "_sizeseq")
	
	C_name2seq = collections.OrderedDict()
	
	for line in files:
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
	
	total_length = 0
	
	for i in C_name2seq:
		seq = C_name2seq[i].rstrip('\n')
		total_length+=len(seq)
		
	
	
	c_length = 0
	
	for i in C_name2seq:
		
		c_length+=len(C_name2seq[i])
		if c_length >= (total_length/2):
			n50 = str(len(C_name2seq[i]))
			break
	return(n50)


def refBlast():
	subprocess.call(["makeblastdb -in " + ref + " -dbtype prot"], shell=True)
	subprocess.call(["blastx -query " + trans + " -db " + ref + " -out transcriptome2refBlastxE30MTS1 -evalue 1e-30 -max_target_seqs 1 -outfmt '6 qseqid sseqid evalue length pident'"], shell=True)
	f = open('transcriptome2refBlastxE30MTS1')
	blastx = f.readlines()
	outblastfile = open('transcriptome2refBlastxE30MTS1_Filtered', 'w')
	for i in blastx:
		length = i.split()[3]
		pident = i.split()[4]
		if int(length) > 100:
			if float(pident) > 50:
				outblastfile.write(i)
	f.close()
	outblastfile.close()

	f = open('transcriptome2refBlastxE30MTS1_Filtered')
	lines = f.readlines()
	refHitCount = 0
	for i in lines:
		refHitCount+=1
	f.close()

	return(refHitCount)




PATH = input(' Give me the full path of the transcriptome you want stats for ')
trans = (PATH)
PATH2 = input('give me the full path of a reference you want to BLAST against ')
ref = (PATH2)

#n50 = n50(trans)

f = open(trans)

X = {}
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

	X[C_name] = C_seq
f.close()


lessthan200 = []
lessthan1kb = []
morethan1kb = []
all_lengths = []


f = open(trans)
lines = f.readlines()
isogroups = []
for i in lines:
	if i.startswith('>'):
		isogroup = i.split('gene=')[1]
		isogroup = isogroup.split('  length')[0]
		isogroups.append(isogroup)

isogroupCount = len(set(isogroups))
f.close()



for i in X:
	all_lengths.append(len(X[i]))
	if len(X[i]) < 200:
		lessthan200.append(i)
	elif len(X[i]) < 1000:
		lessthan1kb.append(i)
	elif len(X[i]) > 1000:
		morethan1kb.append(i)
		
lt200 = int(len(lessthan200) * 100)/len(X)
lt1000 = int(len(lessthan1kb) * 100)/len(X)
mt1000 = int(len(morethan1kb) * 100)/len(X)

refHitCount = refBlast()


print('\n')
print('---------------------------------------------')
print('\n')
print('The total number of sequences is ' + str(len(all_lengths)) + '\n')
print('The number of sequences with a blastx hit to your ref @ e=1e-30 with 50% identity over > 100bp is ' + str(refHitCount) + '\n')
#print('The N50 of ' + trans + ' is ' + n50 + '\n')
print('The total number of isogroups is ' + str(isogroupCount) + '\n')
print('The minumum sequece length is ' +  str(min(all_lengths)) + '\n')
print('The maximum sequece length is ' +  str(max(all_lengths)) + '\n')
print(str(lt200) + '% of sequences are less than 200bp long' + '\n')
print(str(lt1000) + '% of sequences are less than 1000bp long' + '\n')
print(str(mt1000) + '% of sequences are more than 1000bp long' + '\n')
print('\n')
print('---------------------------------------------')







