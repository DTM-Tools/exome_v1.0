from pymongo import MongoClient
import pymongo
import sys
import csv


def getFindingsCollection(ipAddr, srvport):
	client = MongoClient(host=ipAddr, port=srvport)
	db = client.collectionName
	return client.collectionName.tableName

def getFindingsByRunId(findingsCollection, inputF):
    cursor = findingsCollection.find({ "inputFile" : inputF})
    for document in cursor:
        return document

def main(argv):
	inputF = argv[0]
	aFinding = getFindingsByRunId(getFindingsCollection('127.0.0.1', 27017), inputF)

	print ('\n')
	print('DTM-Tools: individual RBC phenotype prediction')
	print('For research use only')
	print()
	print('run_id:', aFinding['runId'])
	print('input file:'+ inputF)
	print ('\n')

	print('**** KEL ****')
	print('Antigens:')
	countHet = 0
	print('K1: ',end='')
	try:
		if (aFinding['findings']['s6a']['zygosity'])== 1 or (aFinding['findings']['s6a']['zygosity'])==2:
			print ('POS')
		elif (aFinding['findings']['s6a']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	print('k(cellano: ',end='')
	try:
		if (aFinding['findings']['s6a']['zygosity'])== 0 or (aFinding['findings']['s6a']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6a']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('Jsa: ',end='')
	try:
		if (aFinding['findings']['s6f']['zygosity'])== 2 or (aFinding['findings']['s6a']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6f']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	print('Jsb: ',end='')
	try:
		if (aFinding['findings']['s6f']['zygosity'])== 0 or (aFinding['findings']['s6f']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6f']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('Ula: ',end='')
	try:
		if (aFinding['findings']['s6g']['zygosity'])== 1 or (aFinding['findings']['s6g']['zygosity'])==2:
			print ('POS')
		elif (aFinding['findings']['s6g']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	print('K12: ',end='')
	try:
		if (aFinding['findings']['s6h']['zygosity'])== 0 or (aFinding['findings']['s6h']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6h']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K11: ',end='')
	try:
		if (aFinding['findings']['s6k']['zygosity'])== 0 or (aFinding['findings']['s6k']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6k']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K17: ',end='')
	try:
		if (aFinding['findings']['s6k']['zygosity'])== 1 or (aFinding['findings']['s6k']['zygosity'])==2:
			print ('POS')
		elif (aFinding['findings']['s6k']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	print('K19: ',end='')
	try:
		if (aFinding['findings']['s6n']['zygosity'])== 0 or (aFinding['findings']['s6n']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6n']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K22: ',end='')
	try:
		if (aFinding['findings']['s6p']['zygosity'])== 0 or (aFinding['findings']['s6p']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6p']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K23: ',end='')
	try:
		if (aFinding['findings']['s6q']['zygosity'])== 1 or (aFinding['findings']['s6q']['zygosity'])==2:
			print ('POS')
		elif (aFinding['findings']['s6q']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	print('K26: ',end='')
	try:
		if (aFinding['findings']['s6s']['zygosity'])== 0 or (aFinding['findings']['s6s']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6s']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K27: ',end='')
	try:
		if (aFinding['findings']['s6t']['zygosity'])== 0 or (aFinding['findings']['s6t']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6t']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K29: ',end='')
	try:
		if (aFinding['findings']['s6v']['zygosity'])== 0 or (aFinding['findings']['s6v']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6v']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K30: ',end='')
	try:
		if (aFinding['findings']['s6w']['zygosity'])== 0 or (aFinding['findings']['s6w']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6w']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K38: ',end='')
	try:
		if (aFinding['findings']['s6x']['zygosity'])== 0 or (aFinding['findings']['s6x']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6x']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K31: ',end='')
	try:
		if (aFinding['findings']['s6x']['zygosity'])== 1 or (aFinding['findings']['s6x']['zygosity'])==2:
			print ('POS')
		elif (aFinding['findings']['s6x']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	print('K32: ',end='')
	try:
		if (aFinding['findings']['s6y']['zygosity'])== 0 or (aFinding['findings']['s6y']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6y']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K33: ',end='')
	try:
		if (aFinding['findings']['s6z']['zygosity'])== 0 or (aFinding['findings']['s6z']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6z']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K34: ',end='')
	try:
		if (aFinding['findings']['s6aa']['zygosity'])== 0 or (aFinding['findings']['s6aa']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6aa']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K36: ',end='')
	try:
		if (aFinding['findings']['s6ab']['zygosity'])== 0 or (aFinding['findings']['s6ab']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6ab']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('K37: ',end='')
	try:
		if (aFinding['findings']['s6ac']['zygosity'])== 0 or (aFinding['findings']['s6ac']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s6ac']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')

	print('Antigens that may require phasing:')
	var1 = 0
	var2 = 0
	print('K14: ',end='')
	try:
		if (aFinding['findings']['s6i']['zygosity'])== 0 or (aFinding['findings']['s6i']['zygosity'])==1 or (aFinding['findings']['s6j']['zygosity'])==0 or (aFinding['findings']['s6j']['zygosity'])==1:
			print ('POS')
			var1 = 1
		elif (aFinding['findings']['s6i']['zygosity'])== 2 or (aFinding['findings']['s6j']['zygosity'])== 2:
			print ('NEG')
		print('K24: ',end='')
		if ((aFinding['findings']['s6j']['zygosity'])== 1 or (aFinding['findings']['s6j']['zygosity'])==2) and ('K24pos' in aFinding['findings']['s6j']['classification']):
			print ('POS')
			var2 = 1
		else:
			print ('NEG')
		if var1 == 0 and var2 ==1:
			print ('Suggest manual review of K14/K24')
	except:
		print('NA')
	print('K18: ',end='')
	try:
		if (aFinding['findings']['s6l']['zygosity'])== 2 or (aFinding['findings']['s6m']['zygosity'])==2:
			print ('NEG')
		else:
			print ('POS')
		if (aFinding['findings']['s6l']['zygosity'])== 1 and (aFinding['findings']['s6m']['zygosity'])==1:
			print ('Suggest reviewing K18 alleles for phasing')
	except:
		print('NA')
	var1=0
	var2=0
	print('K25: ',end='')
	try:
		if (aFinding['findings']['s6r']['zygosity'])== 1 or (aFinding['findings']['s6r']['zygosity'])==2:
			print ('POS')
			var1=1
		elif (aFinding['findings']['s6r']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	print('K28: ',end='')
	try:
		if (aFinding['findings']['s6u']['zygosity'])== 1 or (aFinding['findings']['s6u']['zygosity'])==2:
			print ('POS')
			var2=1
		elif (aFinding['findings']['s6u']['zygosity'])== 0:
			print ('NEG')
		if (var1==1 and var2 ==1) and (aFinding['findings']['s6r']['zygosity']==2 or aFinding['findings']['s6u']['zygosity']==2):
			print ('Suggest reviewing K25/K28 for phasing')
		kpa = 0
		kbp =0
		kpc =0
		if aFinding['findings']['s6o']['zygosity'] == 0:
			print ('Kpc: NEG')
			if aFinding['findings']['s6b']['zygosity'] == 1:
				print ('Kpa: POS')
				print ('Kpb: POS')
			if aFinding['findings']['s6b']['zygosity'] == 2:
				print ('Kpa: POS')
				print ('Kpb: NEG')
			if aFinding['findings']['s6b']['zygosity'] == 0:
				print ('Kpa: NEG')
				print ('Kpb: POS')
		if aFinding['findings']['s6o']['zygosity'] == 1:
			print ('Kpc: POS')
			if aFinding['findings']['s6b']['zygosity'] == 0:
				print ('Kpa: NEG')
				print ('Kpb: POS')
			if aFinding['findings']['s6b']['zygosity'] == 1:
				print ('Kpa: POS')
				print ('Kpb: NEG')
				print ('Suggest phasing review for Kpa/Kpb/Kpc')
			if aFinding['findings']['s6b']['zygosity'] == 2:
				print ('Kpa: POS')
				print ('Kpb: NEG')
				print ('Both Kpa and Kpc appear positive in same haplotype - review manually')
			if aFinding['findings']['s6o']['zygosity'] == 2:
				print ('Kpc: POS')
				if aFinding['findings']['s6b']['zygosity'] == 1 or aFinding['findings']['s6b']['zygosity'] == 2:
					print ('Both Kpa and Kpc appear positive in same haplotype - review Kpa/Kbp/Kpc manually')
				else:
					print('Kpa: NEG')
					print('Kpb: NEG')
	except:
		print('Kpa/Kpb/Kpc NA')
	ListAnt=['s6a','s6b','s6o','s6f','s6g','s6h','s6i','s6j','s6k','s6l','s6m','s6n','s6p','s6r','s6q','s6s','s6u','s6t','s6v','s6w','s6x','s6y','s6z','s6aa','s6ab','s6ac']
	for antigen in ListAnt:
		try:
			if aFinding['findings'][antigen]['zygosity'] == 1:
				countHet=countHet+1
		except:
			pass
	print('Number of heterozygous variants in antigenic sites:',countHet)
	print('\n')
	system = '6'
	idList=[]
	with open('ChromoList.csv') as csvfile:
		linereader = csv.reader(csvfile,delimiter=',')
		for row in linereader:
			if row[3] == system:
				idList.append(row[0])
	with open('ChromoInDelList.csv') as csvfile:
		linereader2 = csv.reader(csvfile,delimiter=',')
		for row2 in linereader2:
			if row2[3] == system:
				idList.append(row2[0])
	with open('File3.csv') as csvfile:
		linereader3 = csv.reader(csvfile,delimiter=',')
		for row3 in linereader3:
			if row3[5] == system:
				idList.append(row3[0])
	countHetWeak=0
	countHomWeak=0
	weakHetList=[]
	weakHomList=[]
	totalWeaks=0
	countHetNull=0
	countHomNull=0
	nullHetList=[]
	nullHomList=[]
	totalNulls=0
	HetMultDoubtWeak=0
	HetMultDoubtWeakList=[]
	HetMultDoubtNull=0
	HetMultDoubtNullList=[]
	for allele in idList:
		try:
			effect = aFinding['findings'][allele]['effect']
			if ('Weak' in effect) and ('part' not in aFinding['findings'][allele]['description']):
				if (aFinding['findings'][allele]['zygosity']==1):
					countHetWeak=countHetWeak+1
					weakHetList.append(aFinding['findings'][allele]['description'])
				else:
					countHomWeak=countHomWeak+1
					weakHomList.append(aFinding['findings'][allele]['description'])
			if ('Null' in effect) and ('part' not in aFinding['findings'][allele]['description']):
				if (aFinding['findings'][allele]['zygosity']==1):
					countHetNull=countHetNull+1
					nullHetList.append(aFinding['findings'][allele]['description'])
				else:
					countHomNull=countHomNull+1
					nullHomList.append(aFinding['findings'][allele]['description'])
		except:
			pass
		try:
			if 'mu' in (aFinding['findings'][allele]):
				try:
					if ('_null_' in aFinding['findings'][allele]['classification']):
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==1):
							countHetNull=countHetNull+1
							nullHetList.append(aFinding['findings'][allele]['description'])
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==0):
							HetMultDoubtNull=HetMultDoubt+1
							HetMultDoubtNullList.append(aFinding['findings'][allele]['description'])
						else:
							countHomNull=countHomNull+1
							nullHomList.append(aFinding['findings'][allele]['description'])
				except:
					pass
		except:
			pass
		try:
			if ('mu' in aFinding['findings'][allele]):
				try:
					if ('_weak_' in aFinding['findings'][allele]['classification']):
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==1):
							countHetWeak=countHetWeak+1
							weakHetList.append(aFinding['findings'][allele]['description'])
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==0):
							HetMultDoubtWeak=HetMultDoubt+1
							HetMultDoubtWeakList.append(aFinding['findings'][allele]['description'])
						else:
							countHomWeak=countHomWeak+1
							weakHomList.append(aFinding['findings'][allele]['description'])
				except:
					pass
		except:
			pass
	totalNulls=countHetNull + countHomNull + HetMultDoubtWeak
	totalWeaks=countHetWeak + countHomWeak + HetMultDoubtNull
	print('Total Number of weak alleles found: ',totalWeaks)
	if totalWeaks > 0:
		print('\t'+'Number of heterozygous weak alleles found: ',countHetWeak)
		if countHetWeak >0:
			print('\t'+str(weakHetList))
		print('\t'+'Number of homozygous weak alleles found: ',countHomWeak)
		if countHomWeak >0:
			print('\t'+str(weakHomList))
		print('\t'+'Number of homozygous weak alleles defined by more than one SNV/indel possible (certainty of zero) ',HetMultDoubtWeak)
		if HetMultDoubtWeak >0:
			print('\t'+str(HetMultDoubtWeakList))

	print('\n')
	print('Total Number of null alleles found: ',totalNulls)
	if totalNulls > 0:
		print('\t'+'Total number of heterozygous null alleles found: ',countHetNull)
		if countHetNull >0:
			print('\t'+str(nullHetList))
		print('\t'+'Number of homozygous null alleles found: ',countHomNull)
		if countHomNull >0:
			print('\t'+str(nullHomList))
		print('\t'+'Number of heterozygous null alleles defined by more than one SNV/indel possible (certainty of zero) ',HetMultDoubtNull)
		if HetMultDoubtNull >0:
			print('\t'+str(HetMultDoubtNullList))
	print ('\n')

	print ('Filtered alleles: ', end = '')
	filtered = 0
	for id in idList:
		try:
			if aFinding['findings'][id]['determination'] == 'filtered':
				print (aFinding['findings'][id]['description'],'; ',end='')
				filtered = filtered + 1
		except:
			pass
	if filtered == 0:
		print ('None')
	print ('\n')

	print ('Regions without aligned reads: ', end = '')
	missing = 0
	for id in idList:
		try:
			if aFinding['findings'][id]['determination'] == 'MISSING':
				print (aFinding['findings'][id]['description'],'; ',end='')
				missing = missing + 1
		except:
			pass
	if missing == 0:
		print ('None')
	print ('\n')

	print ('Additional variants in KEL exons:')
	novel = (aFinding['additionalFindings'])
	for sample in novel.keys():
		if novel[sample]['gene'] == 'KEL':
			if novel[sample]['location'] == 'exonic':
				print (novel[sample]['identifier'], ':', novel[sample]['varType'])
	print('\n')

	print('**** DUFFY ****')
	print('Antigens:')
	countHet = 0
	print('Fya: ',end='')
	try:
		if (aFinding['findings']['s8a']['zygosity'])== 0 or (aFinding['findings']['s8a']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s8a']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('Fyb: ',end='')
	try:
		if (aFinding['findings']['s8a']['zygosity'])== 2 or (aFinding['findings']['s8a']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s8a']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	try:
		if (aFinding['findings']['s8a']['zygosity'])== 1:
			countHet=1
	except:
		pass
	print('Number of heterozygous antigenic calls:',countHet)
	print('\n')
	system = '8'
	idList=[]
	with open('ChromoList.csv') as csvfile:
		linereader = csv.reader(csvfile,delimiter=',')
		for row in linereader:
			if row[3] == system:
				idList.append(row[0])
	with open('ChromoInDelList.csv') as csvfile:
		linereader2 = csv.reader(csvfile,delimiter=',')
		for row2 in linereader2:
			if row2[3] == system:
				idList.append(row2[0])
	with open('File3.csv') as csvfile:
		linereader3 = csv.reader(csvfile,delimiter=',')
		for row3 in linereader3:
			if row3[5] == system:
				idList.append(row3[0])
	countHetWeak=0
	countHomWeak=0
	weakHetList=[]
	weakHomList=[]
	totalWeaks=0
	countHetNull=0
	countHomNull=0
	nullHetList=[]
	nullHomList=[]
	totalNulls=0
	HetMultDoubtWeak=0
	HetMultDoubtWeakList=[]
	HetMultDoubtNull=0
	HetMultDoubtNullList=[]
	for allele in idList:
		try:
			effect = aFinding['findings'][allele]['effect']
			if ('Weak' in effect) and ('part' not in aFinding['findings'][allele]['description']):
				if (aFinding['findings'][allele]['zygosity']==1):
					countHetWeak=countHetWeak+1
					weakHetList.append(aFinding['findings'][allele]['description'])
				else:
					countHomWeak=countHomWeak+1
					weakHomList.append(aFinding['findings'][allele]['description'])
			if ('Null' in effect) and ('part' not in aFinding['findings'][allele]['description']):
				if (aFinding['findings'][allele]['zygosity']==1):
					countHetNull=countHetNull+1
					nullHetList.append(aFinding['findings'][allele]['description'])
				else:
					countHomNull=countHomNull+1
					nullHomList.append(aFinding['findings'][allele]['description'])
		except:
			pass
		try:
			if ('mu' in aFinding['findings'][allele]):
				try:
					if ('_null_' in aFinding['findings'][allele]['classification']):
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==1):
							countHetNull=countHetNull+1
							nullHetList.append(aFinding['findings'][allele]['description'])
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==0):
							HetMultDoubtNull=HetMultDoubt+1
							HetMultDoubtNullList.append(aFinding['findings'][allele]['description'])
						else:
							countHomNull=countHomNull+1
							nullHomList.append(aFinding['findings'][allele]['description'])
				except:
					pass
		except:
			pass
		try:
			if ('mu' in aFinding['findings'][allele]):
				try:
					if ('_weak_' in aFinding['findings'][allele]['classification']):
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==1):
							countHetWeak=countHetWeak+1
							weakHetList.append(aFinding['findings'][allele]['description'])
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==0):
							HetMultDoubtWeak=HetMultDoubt+1
							HetMultDoubtWeakList.append(aFinding['findings'][allele]['description'])
						else:
							countHomWeak=countHomWeak+1
							weakHomList.append(aFinding['findings'][allele]['description'])
				except:
					pass
		except:
			pass
	totalNulls=countHetNull + countHomNull + HetMultDoubtWeak
	totalWeaks=countHetWeak + countHomWeak + HetMultDoubtNull
	print('Total Number of weak alleles found: ',totalWeaks)
	if totalWeaks > 0:
		print('\t'+'Number of heterozygous weak alleles found: ',countHetWeak)
		if countHetWeak >0:
			print('\t'+str(weakHetList))
		print('\t'+'Number of homozygous weak alleles found: ',countHomWeak)
		if countHomWeak >0:
			print('\t'+str(weakHomList))
		print('\t'+'Number of homozygous weak alleles defined by more than one SNV/indel possible (certainty of zero) ',HetMultDoubtWeak)
		if HetMultDoubtWeak >0:
			print('\t'+str(HetMultDoubtWeakList))

	print('\n')
	print('Total Number of null alleles found: ',totalNulls)
	if totalNulls > 0:
		print('\t'+'Total number of heterozygous null alleles found: ',countHetNull)
		if countHetNull >0:
			print('\t'+str(nullHetList))
		print('\t'+'Number of homozygous null alleles found: ',countHomNull)
		if countHomNull >0:
			print('\t'+str(nullHomList))
		print('\t'+'Number of heterozygous null alleles defined by more than one SNV/indel possible (certainty of zero) ',HetMultDoubtNull)
		if HetMultDoubtNull >0:
			print('\t'+str(HetMultDoubtNullList))
	print ('\n')

	print ('Filtered alleles: ', end = '')
	filtered = 0
	for id in idList:
		try:
			if aFinding['findings'][id]['determination'] == 'filtered':
				print (aFinding['findings'][id]['description'],'; ',end='')
				filtered = filtered + 1
		except:
			pass
	if filtered == 0:
		print ('None')
	print ('\n')

	print ('Regions without aligned reads: ', end = '')
	missing = 0
	for id in idList:
		try:
			if aFinding['findings'][id]['determination'] == 'MISSING':
				print (aFinding['findings'][id]['description'],'; ',end='')
				missing = missing + 1
		except:
			pass
	if missing == 0:
		print ('None')
	print ('\n')

	print ('Additional variants in ACKR1 exons:')
	novel = (aFinding['additionalFindings'])
	for sample in novel.keys():
		if novel[sample]['gene'] == 'ACKR1':
			if novel[sample]['location'] == 'exonic':
				print (novel[sample]['identifier'], ':', novel[sample]['varType'])
	print('\n')

	print('**** KIDD ****')
	print('Antigens:')
	countHet = 0
	print('Jka: ',end='')
	try:
		if (aFinding['findings']['s9a']['zygosity'])== 0 or (aFinding['findings']['s9a']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s9a']['zygosity'])== 2:
			print ('NEG')
	except:
		print('NA')
	print('Jkb: ',end='')
	try:
		if (aFinding['findings']['s9a']['zygosity'])== 2 or (aFinding['findings']['s9a']['zygosity'])==1:
			print ('POS')
		elif (aFinding['findings']['s9a']['zygosity'])== 0:
			print ('NEG')
	except:
		print('NA')
	try:
		if (aFinding['findings']['s9a']['zygosity'])== 1:
			countHet=1
	except:
		pass
	print('Number of heterozygous antigenic calls:',countHet)
	print('\n')
	system = '9'
	idList=[]
	with open('ChromoList.csv') as csvfile:
		linereader = csv.reader(csvfile,delimiter=',')
		for row in linereader:
			if row[3] == system:
				idList.append(row[0])
	with open('ChromoInDelList.csv') as csvfile:
		linereader2 = csv.reader(csvfile,delimiter=',')
		for row2 in linereader2:
			if row2[3] == system:
				idList.append(row2[0])
	with open('File3.csv') as csvfile:
		linereader3 = csv.reader(csvfile,delimiter=',')
		for row3 in linereader3:
			if row3[5] == system:
				idList.append(row3[0])
	countHetWeak=0
	countHomWeak=0
	weakHetList=[]
	weakHomList=[]
	totalWeaks=0
	countHetNull=0
	countHomNull=0
	nullHetList=[]
	nullHomList=[]
	totalNulls=0
	HetMultDoubtWeak=0
	HetMultDoubtWeakList=[]
	HetMultDoubtNull=0
	HetMultDoubtNullList=[]
	for allele in idList:
		try:
			effect = aFinding['findings'][allele]['effect']
			if ('Weak' in effect) and ('part' not in aFinding['findings'][allele]['description']):
				if (aFinding['findings'][allele]['zygosity']==1):
					countHetWeak=countHetWeak+1
					weakHetList.append(aFinding['findings'][allele]['description'])
				else:
					countHomWeak=countHomWeak+1
					weakHomList.append(aFinding['findings'][allele]['description'])
			if ('Null' in effect) and ('part' not in aFinding['findings'][allele]['description']):
				if (aFinding['findings'][allele]['zygosity']==1):
					countHetNull=countHetNull+1
					nullHetList.append(aFinding['findings'][allele]['description'])
				else:
					countHomNull=countHomNull+1
					nullHomList.append(aFinding['findings'][allele]['description'])
		except:
			pass
		try:
			if ('mu' in aFinding['findings'][allele]):
				try:
					if ('_null_' in aFinding['findings'][allele]['classification']):
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==1):
							countHetNull=countHetNull+1
							nullHetList.append(aFinding['findings'][allele]['description'])
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==0):
							HetMultDoubtNull=HetMultDoubt+1
							HetMultDoubtNullList.append(aFinding['findings'][allele]['description'])
						else:
							countHomNull=countHomNull+1
							nullHomList.append(aFinding['findings'][allele]['description'])
				except:
					pass
		except:
			pass
		try:
			if ('mu' in aFinding['findings'][allele]):
				try:
					if ('_weak_' in aFinding['findings'][allele]['classification']):
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==1):
							countHetWeak=countHetWeak+1
							weakHetList.append(aFinding['findings'][allele]['description'])
						if (aFinding['findings'][allele]['zygosity']==1) and (aFinding['findings'][allele]['certainty']==0):
							HetMultDoubtWeak=HetMultDoubt+1
							HetMultDoubtWeakList.append(aFinding['findings'][allele]['description'])
						else:
							countHomWeak=countHomWeak+1
							weakHomList.append(aFinding['findings'][allele]['description'])
				except:
					pass
		except:
			pass
	totalNulls=countHetNull + countHomNull + HetMultDoubtWeak
	totalWeaks=countHetWeak + countHomWeak + HetMultDoubtNull
	print('Total Number of weak alleles found: ',totalWeaks)
	if totalWeaks > 0:
		print('\t'+'Number of heterozygous weak alleles found: ',countHetWeak)
		if countHetWeak >0:
			print('\t'+str(weakHetList))
		print('\t'+'Number of homozygous weak alleles found: ',countHomWeak)
		if countHomWeak >0:
			print('\t'+str(weakHomList))
		print('\t'+'Number of homozygous weak alleles defined by more than one SNV/indel possible (certainty of zero) ',HetMultDoubtWeak)
		if HetMultDoubtWeak >0:
			print('\t'+str(HetMultDoubtWeakList))

	print('\n')
	print('Total Number of null alleles found: ',totalNulls)
	if totalNulls > 0:
		print('\t'+'Total number of heterozygous null alleles found: ',countHetNull)
		if countHetNull >0:
			print('\t'+str(nullHetList))
		print('\t'+'Number of homozygous null alleles found: ',countHomNull)
		if countHomNull >0:
			print('\t'+str(nullHomList))
		print('\t'+'Number of heterozygous null alleles defined by more than one SNV/indel possible (certainty of zero) ',HetMultDoubtNull)
		if HetMultDoubtNull >0:
			print('\t'+str(HetMultDoubtNullList))
	print ('\n')

	print ('Filtered alleles: ', end = '')
	filtered = 0
	for id in idList:
		try:
			if aFinding['findings'][id]['determination'] == 'filtered':
				print (aFinding['findings'][id]['description'],'; ',end='')
				filtered = filtered + 1
		except:
			pass
	if filtered == 0:
		print ('None')
	print ('\n')

	print ('Regions without aligned reads: ', end = '')
	missing = 0
	for id in idList:
		try:
			if aFinding['findings'][id]['determination'] == 'MISSING':
				print (aFinding['findings'][id]['description'],'; ',end='')
				missing = missing + 1
		except:
			pass
	if missing == 0:
		print ('None')
	print ('\n')

	print ('Additional variants in SLC14A1 exons:')
	novel = (aFinding['additionalFindings'])
	for sample in novel.keys():
		if novel[sample]['gene'] == 'SLC14A1':
			if novel[sample]['location'] == 'exonic':
				print (novel[sample]['identifier'], ':', novel[sample]['varType'])
	print('\n')

if __name__ == "__main__":
   main(sys.argv[1:])
