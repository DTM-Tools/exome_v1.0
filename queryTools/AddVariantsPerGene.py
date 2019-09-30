from pymongo import MongoClient
import pymongo
import sys
import csv
import statistics

def getNovelList(ipAddr, srvport):
    client = MongoClient(host=ipAddr, port=srvport)
    alleles = {}
    novel = client.collectionName.tableName
    foundNovel = novel.find()

    for person in foundNovel:
        alleles.update(person['additionalFindings'])

    client.close()
    return alleles

alleles = getNovelList('127.0.0.1', 27017)
countExonic = 0
nonsyn = 0
synon = 0
nonframesubs = 0
nonframeIns=0
nonframeDel=0
frameDel =0
frameIns=0
frameSubs =0
newStop=0
lossStop=0
unknowns=0
idList = []
varTypeList = []
for id in alleles.keys():
    if id not in idList:
        if alleles[id]['location']=="exonic":
            countExonic = countExonic + 1
            idList.append(id)
            if alleles[id]['varType'] not in varTypeList:
                varTypeList.append(alleles[id]['varType'])
            if alleles[id]['varType']=='nonsynonymous SNV':
                nonsyn=nonsyn+1
            elif alleles[id]['varType']=='synonymous SNV':
                synon=synon+1
            elif alleles[id]['varType']=='nonframeshift substitution':
                nonframesubs = nonframesubs +1
            elif alleles[id]['varType']=='nonframeshift insertion':
                nonframeIns = nonframeIns +1
            elif alleles[id]['varType']=='nonframeshift deletion':
                nonframeDel = nonframeDel +1
            elif alleles[id]['varType']=='frameshift deletion':
                frameDel = frameDel +1
            elif alleles[id]['varType']=='frameshift insertion':
                frameIns = frameIns +1
            elif alleles[id]['varType']=='frameshift substitution':
                frameSubs = frameSubs +1
            elif alleles[id]['varType']=='stopgain':
                newStop = newStop +1
            elif alleles[id]['varType']=='stoploss':
                lossStop = lossStop+1
            elif alleles[id]['varType']=='unknown':
                unknowns = unknowns +1

print (varTypeList)
print ('Total exonic changes:', + countExonic)
print ('Total nonsynonymous SNV changes:', + nonsyn)
print ('Total synonymous SNV changes:', + synon)
print ('Total nonframeshift substitution:', + nonframesubs)
print ('Total nonframeshift insertion:', + nonframeIns)
print ('Total nonframeshift deletion:', + nonframeDel)
print ('Total frameshift deletion:', + frameDel)
print ('Total frameshift insertion:', + frameIns)
print ('Total frameshift substitution:', + frameSubs)
print ('Total stop gains:', + newStop)
print ('Total stop loss:', + lossStop)
print ('Total exonic unknowns:', + unknowns)

idListIntronic=[]
countIntronic=0
for id in alleles.keys():
    if id not in idListIntronic:
        if alleles[id]['location']=="intronic":
            countIntronic = countIntronic + 1
            idListIntronic.append(id)
print ('Total intronic changes:',+countIntronic)

idListUTR3=[]
countUTR3=0
for id in alleles.keys():
    if id not in idListUTR3:
        if alleles[id]['location']=="UTR3":
            countUTR3 = countUTR3 + 1
            idListUTR3.append(id)
print ('Total UTR3 changes:',+countUTR3)

idListUTR5=[]
countUTR5=0
for id in alleles.keys():
    if id not in idListUTR5:
        if alleles[id]['location']=="UTR5":
            countUTR5 = countUTR5 + 1
            idListUTR5.append(id)
print ('Total UTR5 changes:',+countUTR5)

geneList=[]
countGenes=0
for id in alleles.keys():
    gene = alleles[id]['gene']
    geneCut = gene.split('(')[0]
    if geneCut not in geneList:
        geneList.append(geneCut)
        countGenes=countGenes+1
print ('Total genes:',+countGenes)

for element in geneList:
    print()
    print('Stats per gene:')
    print('*** '+str(element),'***')
    countExonic = 0
    nonsyn = 0
    synon = 0
    nonframesubs = 0
    nonframeIns=0
    nonframeDel=0
    frameDel =0
    frameIns=0
    frameSubs =0
    newStop=0
    lossStop=0
    unknowns=0
    geneSPDIList=[]
    for id in alleles.keys():
        if alleles[id]['gene']==element:
            if id not in geneSPDIList:
                geneSPDIList.append(id)
                if alleles[id]['location']=="exonic":
                    countExonic = countExonic + 1
                    if alleles[id]['varType']=='nonsynonymous SNV':
                        nonsyn=nonsyn+1
                    elif alleles[id]['varType']=='synonymous SNV':
                        synon=synon+1
                    elif alleles[id]['varType']=='nonframeshift substitution':
                        nonframesubs = nonframesubs +1
                    elif alleles[id]['varType']=='nonframeshift insertion':
                        nonframeIns = nonframeIns +1
                    elif alleles[id]['varType']=='nonframeshift deletion':
                        nonframeDel = nonframeDel +1
                    elif alleles[id]['varType']=='frameshift deletion':
                        frameDel = frameDel +1
                    elif alleles[id]['varType']=='frameshift insertion':
                        frameIns = frameIns +1
                    elif alleles[id]['varType']=='frameshift substitution':
                        frameSubs = frameSubs +1
                    elif alleles[id]['varType']=='stopgain':
                        newStop = newStop +1
                    elif alleles[id]['varType']=='stoploss':
                        lossStop = lossStop+1
                    elif alleles[id]['varType']=='unknown':
                        unknowns = unknowns +1

    print ('Total exonic changes:', + countExonic)
    print ('Total nonsynonymous SNV changes:', + nonsyn)
    print ('Total synonymous SNV changes:', + synon)
    print ('Total nonframeshift substitution:', + nonframesubs)
    print ('Total nonframeshift insertion:', + nonframeIns)
    print ('Total nonframeshift deletion:', + nonframeDel)
    print ('Total frameshift deletion:', + frameDel)
    print ('Total frameshift insertion:', + frameIns)
    print ('Total frameshift substitution:', + frameSubs)
    print ('Total stop gains:', + newStop)
    print ('Total stop loss:', + lossStop)
    print ('Total exonic unknowns:', + unknowns)

    idListIntronic=[]
    countIntronic=0
    for id in alleles.keys():
        if id not in idListIntronic:
            if alleles[id]['location']=="intronic":
                countIntronic = countIntronic + 1
                idListIntronic.append(id)
    print ('Total intronic changes:',+countIntronic)

    idListUTR3=[]
    countUTR3=0
    for id in alleles.keys():
        if id not in idListUTR3:
            if alleles[id]['location']=="UTR3":
                countUTR3 = countUTR3 + 1
                idListUTR3.append(id)
    print ('Total UTR3 changes:',+countUTR3)

    idListUTR5=[]
    countUTR5=0
    for id in alleles.keys():
        if id not in idListUTR5:
            if alleles[id]['location']=="UTR5":
                countUTR5 = countUTR5 + 1
                idListUTR5.append(id)
    print ('Total UTR5 changes:',+countUTR5)
