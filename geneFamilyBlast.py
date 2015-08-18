hitLists = []
f = open('CON_unique_fastaSelfBlastNoE')
blast = f.readlines()
for i in blast:
	q = i.split()[0]
	d = i.split()[1]
	l = i.split()[2]
	p = i.split()[3]
	e = i.split()[4]
	if q != d:
		if int(l) > 150:
			if float(p) > 40:
				hitPair = [q,d]
				hitLists.append(hitPair)

merged = []
outFile = open('clusters', 'w')
for lists in [set(d) for d in hitLists]:
	for m in merged:
		if not m.isdisjoint(lists):
			m.update(lists)
			break
	else:
		merged.append(lists)
for i in merged:
	for s in i:
		outFile.write(s.rstrip('\n') + '\t')
	outFile.write('\n')
















