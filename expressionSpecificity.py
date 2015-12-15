# Script to get stats on expression per cluster
# Presence/absence is at a threshold of TPM > 2
# 



def percentFloat(proportion,total):
	percentage = (float(proportion) * 100)/float(total)
	return(percentage)


def exp(clust):
	myclusts = clust.split()
	
	# Stores percentages of tissues transcripts in cluster are expressed in
	CONclusterOnPercentages = []
	CONclusterOffPercentages = []
	PLEclusterOnPercentages = []
	PLEclusterOffPercentages = []
	
	# Stores names of transcripts which are expressed in these tissues
	CONff = []
	CONleaf = []
	CONmf = []
	CONpet = []
	CONroot = []
	CONvb = []
	# Stores expression levels of transcripts which are expressed in these tissues
	CONffExpression = []
	CONleafExpression = []
	CONmfExpression = []
	CONpetExpression = []
	CONrootExpression = []
	CONvbExpression = []
	# Stores names of transcripts which are expressed in these tissues
	PLEff = []
	PLEleaf = []
	PLEmf = []
	PLEpet = []
	PLEroot = []
	PLEvb = []
	# Stores expression levels of transcripts which are expressed in these tissues
	PLEffExpression = []
	PLEleafExpression = []
	PLEmfExpression = []
	PLEpetExpression = []
	PLErootExpression = []
	PLEvbExpression = []
	CONperCopyExpression = []
	PLEperCopyExpression = []

	for cl in myclusts:
		
		CONon = [] 
		CONoff = []
		CONexpression = []
		CONtotal = []
		
		
		PLEon = []
		PLEoff = []
		PLEexpression = []
		PLEtotal = []
		
		if cl.startswith('CON'):
			cl = cl.lstrip('CON')
			if float((CONtpmDict[cl])[0]) > 2:
				CONff.append(cl)
				CONffExpression.append((CONtpmDict[cl])[0])
				CONon.append(cl)
				CONexpression.append((CONtpmDict[cl])[0])
				CONtotal.append(cl)
			elif float((CONtpmDict[cl])[0]) < 2:
				CONoff.append(cl)
				CONtotal.append(cl)
			
			if float((CONtpmDict[cl])[1]) > 2:
				CONleaf.append(cl)
				CONleafExpression.append((CONtpmDict[cl])[1])
				CONon.append(cl)
				CONexpression.append((CONtpmDict[cl])[1])
				CONtotal.append(cl)
			elif float((CONtpmDict[cl])[1]) < 2:
				CONoff.append(cl)
				CONtotal.append(cl)
			
			if float((CONtpmDict[cl])[2]) > 2:
				CONmf.append(cl)
				CONmfExpression.append((CONtpmDict[cl])[2])
				CONon.append(cl)
				CONexpression.append((CONtpmDict[cl])[2])
				CONtotal.append(cl)
			elif float((CONtpmDict[cl])[2]) < 2:
				CONoff.append(cl)
				CONtotal.append(cl)
			
			if float((CONtpmDict[cl])[3]) > 2:
				CONpet.append(cl)
				CONpetExpression.append((CONtpmDict[cl])[3])
				CONon.append(cl)
				CONexpression.append((CONtpmDict[cl])[3])
				CONtotal.append(cl)
			elif float((CONtpmDict[cl])[3]) < 2:
				CONoff.append(cl)
				CONtotal.append(cl)
			
			if float((CONtpmDict[cl])[4]) > 2:
				CONroot.append(cl)
				CONrootExpression.append((CONtpmDict[cl])[4])
				CONon.append(cl)
				CONexpression.append((CONtpmDict[cl])[4])
				CONtotal.append(cl)
			elif float((CONtpmDict[cl])[4]) < 2:
				CONoff.append(cl)
				CONtotal.append(cl)
			
			if float((CONtpmDict[cl])[5]) > 2:
				CONvb.append(cl)
				CONvbExpression.append((CONtpmDict[cl])[5])
				CONon.append(cl)
				CONexpression.append((CONtpmDict[cl])[5])
				CONtotal.append(cl)
			elif float((CONtpmDict[cl])[5]) < 2:
				CONoff.append(cl)
				CONtotal.append(cl)

			CONmeanExpression = ''
			CONminExpression = ''
			CONmaxExpression = ''
			if CONexpression:
				CONmeanExpression = (sum(map(float,CONexpression))/int(len(CONexpression)))
				CONminExpression = (min(map(float,CONexpression)))
				CONmaxExpression = (max(map(float,CONexpression)))
			else:
				CONmeanExpression = 'NA'
				CONminExpression = 'NA'
				CONmaxExpression = 'NA'
			CONpercentOn = percentFloat((len(CONon)),(len(CONtotal)))
			CONclusterOnPercentages.append(CONpercentOn)
			CONoffCount = (len(CONoff))
			CONonCount = (len(CONon))
			
			# Store per-copy expression information in a string, append to a list (to be returned at the end of function)

			CONperCopy = (cl + '\t' + str(CONmeanExpression) + '\t' + str(CONminExpression) + '\t' + str(CONmaxExpression) + '\t' + str(CONpercentOn))
			CONperCopyExpression.append(CONperCopy)


		if cl.startswith('PLE'):
			cl = cl.lstrip('PLE')
			if float((PLEtpmDict[cl])[0]) > 2:
				PLEff.append(cl)
				PLEffExpression.append((PLEtpmDict[cl])[0])
				PLEon.append(cl)
				PLEexpression.append((PLEtpmDict[cl])[0])
				PLEtotal.append(cl)
			elif float((PLEtpmDict[cl])[0]) < 2:
				PLEoff.append(cl)
				PLEtotal.append(cl)

			if float((PLEtpmDict[cl])[1]) > 2:
				PLEleaf.append(cl)
				PLEleafExpression.append((PLEtpmDict[cl])[1])
				PLEon.append(cl)
				PLEexpression.append((PLEtpmDict[cl])[1])
				PLEtotal.append(cl)
			elif float((PLEtpmDict[cl])[1]) < 2:
				PLEoff.append(cl)
				PLEtotal.append(cl)

			if float((PLEtpmDict[cl])[2]) > 2:
				PLEmf.append(cl)
				PLEmfExpression.append((PLEtpmDict[cl])[2])
				PLEon.append(cl)
				PLEexpression.append((PLEtpmDict[cl])[2])
				PLEtotal.append(cl)
			elif float((PLEtpmDict[cl])[2]) < 2:
				PLEoff.append(cl)
				PLEtotal.append(cl)

			if float((PLEtpmDict[cl])[3]) > 2:
				PLEpet.append(cl)
				PLEpetExpression.append((PLEtpmDict[cl])[3])
				PLEon.append(cl)
				PLEexpression.append((PLEtpmDict[cl])[3])
				PLEtotal.append(cl)
			elif float((PLEtpmDict[cl])[3]) < 2:
				PLEoff.append(cl)
				PLEtotal.append(cl)

			if float((PLEtpmDict[cl])[4]) > 2:
				PLEroot.append(cl)
				PLErootExpression.append((PLEtpmDict[cl])[4])
				PLEon.append(cl)
				PLEexpression.append((PLEtpmDict[cl])[4])
				PLEtotal.append(cl)
			elif float((PLEtpmDict[cl])[4]) < 2:
				PLEoff.append(cl)
				PLEtotal.append(cl)

			if float((PLEtpmDict[cl])[5]) > 2:
				PLEvb.append(cl)
				PLEvbExpression.append((PLEtpmDict[cl])[5])
				PLEon.append(cl)
				PLEexpression.append((PLEtpmDict[cl])[5])
				PLEtotal.append(cl)
			elif float((PLEtpmDict[cl])[5]) < 2:
				PLEoff.append(cl)
				PLEtotal.append(cl)

			PLEmeanExpression = ''
			PLEminExpression = ''
			PLEmaxExpression = ''
			if PLEexpression:
				PLEmeanExpression = (sum(map(float,PLEexpression))/int(len(PLEexpression)))
				PLEminExpression = (min(map(float,PLEexpression)))
				PLEmaxExpression = (max(map(float,PLEexpression)))
			else:
				PLEmeanExpression = 'NA'
				PLEminExpression = 'NA'
				PLEmaxExpression = 'NA'
			PLEpercentOn = percentFloat((len(PLEon)),(len(PLEtotal)))
			PLEclusterOnPercentages.append(PLEpercentOn)
			PLEoffCount = (len(PLEoff))
			PLEonCount = (len(PLEon))

			PLEperCopy = (cl + '\t' + str(PLEmeanExpression) + '\t' + str(PLEminExpression) + '\t' + str(PLEmaxExpression) + '\t' + str(PLEpercentOn))
			PLEperCopyExpression.append(PLEperCopy)


	CONffMeanExpression = ''
	CONleafMeanExpression = ''
	CONmfMeanExpression = ''
	CONpetMeanExpression = ''
	CONrootMeanExpression = ''
	CONvbMeanExpression = ''
	PLEffMeanExpression = ''
	PLEleafMeanExpression = ''
	PLEmfMeanExpression = ''
	PLEpetMeanExpression = ''
	PLErootMeanExpression = ''
	PLEvbMeanExpression = ''
	

	if CONffExpression:
		CONffMeanExpression = (sum(map(float,CONffExpression))/int(len(CONffExpression)))
	else:
		CONffMeanExpression = 'NA'
	if CONleafExpression:
		CONleafMeanExpression = (sum(map(float,CONleafExpression))/int(len(CONleafExpression)))
	else:
		CONleafMeanExpression = 'NA'
	if CONmfExpression:
		CONmfMeanExpression = (sum(map(float,CONmfExpression))/int(len(CONmfExpression)))
	else:
		CONmfMeanExpression = 'NA'
	if CONpetExpression:
		CONpetMeanExpression = (sum(map(float,CONpetExpression))/int(len(CONpetExpression)))
	else:
		CONpetMeanExpression = 'NA'
	if CONrootExpression:
		CONrootMeanExpression = (sum(map(float,CONrootExpression))/int(len(CONrootExpression)))
	else:
		CONrootMeanExpression = 'NA'
	if CONvbExpression:
		CONvbMeanExpression = (sum(map(float,CONvbExpression))/int(len(CONvbExpression)))
	else:
		CONvbMeanExpression = 'NA'


	if PLEffExpression:
		PLEffMeanExpression = (sum(map(float,PLEffExpression))/int(len(PLEffExpression)))
	else:
		PLEffMeanExpression = 'NA'
	if PLEleafExpression:
		PLEleafMeanExpression = (sum(map(float,PLEleafExpression))/int(len(PLEleafExpression)))
	else:
		PLEleafMeanExpression = 'NA'
	if PLEmfExpression:
		PLEmfMeanExpression = (sum(map(float,PLEmfExpression))/int(len(PLEmfExpression)))
	else:
		PLEmfMeanExpression = 'NA'
	if PLEpetExpression:
		PLEpetMeanExpression = (sum(map(float,PLEpetExpression))/int(len(PLEpetExpression)))
	else:
		PLEpetMeanExpression = 'NA'
	if PLErootExpression:
		PLErootMeanExpression = (sum(map(float,PLErootExpression))/int(len(PLErootExpression)))
	else:
		PLErootMeanExpression = 'NA'
	if PLEvbExpression:
		PLEvbMeanExpression = (sum(map(float,PLEvbExpression))/int(len(PLEvbExpression)))
	else:
		PLEvbMeanExpression = 'NA'

	CONpresenceAbsence = ('CONpresenceAbsence' + '\t' + str(len(CONff)) + '\t' + str(len(CONleaf)) + '\t' + str(len(CONmf)) + '\t' + str(len(CONpet)) + '\t' + str(len(CONroot)) + '\t' + str(len(CONvb)))
	CONmeanExpression = ('CONMeanExpression' + '\t' + str(CONffMeanExpression) + '\t' + str(CONleafMeanExpression) + '\t' + str(CONmfMeanExpression) + '\t' + str(CONpetMeanExpression) + '\t' + str(CONrootMeanExpression) + '\t' + str(CONvbMeanExpression))
	PLEpresenceAbsence = ('PLEpresenceAbsence' + '\t' + str(len(PLEff)) + '\t' + str(len(PLEleaf)) + '\t' + str(len(PLEmf)) + '\t' + str(len(PLEpet)) + '\t' + str(len(PLEroot)) + '\t' + str(len(PLEvb)))
	PLEmeanExpression = ('PLEMeanExpression' + '\t' + str(PLEffMeanExpression) + '\t' + str(PLEleafMeanExpression) + '\t' + str(PLEmfMeanExpression) + '\t' + str(PLEpetMeanExpression) + '\t' + str(PLErootMeanExpression) + '\t' + str(PLEvbMeanExpression))

	return CONpresenceAbsence, CONmeanExpression, PLEpresenceAbsence, PLEmeanExpression, CONclusterOnPercentages, PLEclusterOnPercentages, CONperCopyExpression, PLEperCopyExpression



CONtpmDict = {}
c = open('CONmeanTPMs')
c = c.readlines()
for i in c:
	i = i.lstrip()
	if not i.startswith('Trimmed'):
		ref = i.split()[0]
		ref = ref.split('_')[0]
		tissues = i.split()[1:7]
		CONtpmDict[ref] = (tissues)

PLEtpmDict = {}
p = open('PLEmeanTPMs')
p = p.readlines()
for j in p:
	j = j.lstrip()
	if not j.startswith('Trimmed'):
		ref = j.split()[0]
		ref = ref.split('_')[0]
		tissues = j.split()[1:7]
		PLEtpmDict[ref] = (tissues)



f = open('clusters')
clusts = f.readlines()
for clust in clusts:
	CONpa, CONme, PLEpa, PLEme, CONclustOnPercent, PLEclustOnPercent, CONpce, PLEpce = exp(clust)

	
	








