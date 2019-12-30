from pymongo import MongoClient
import pymongo
import sys
import csv

def getFindingsCollection(ipAddr, srvport):
	client = MongoClient(host=ipAddr, port=srvport)
	db = client.rylanFindings
	return client.collectionName.tableName

def getCountOfMissing(findingscollection, allele):
	cursor = findingscollection.count({"findings." + allele + ".determination" : "MISSING"})
	print ("The number of individuals with this locus MISSING in input file is " + str(cursor))
	cursor = findingscollection.count({"findings." + allele + ".determination" : "filtered"})
	print ("The number of individuals with this locus FILTERED is " + str(cursor))

def getCountOfAlleles(findingscollection, allele, zygosity, description):
	cursor = findingscollection.count({"findings." + allele + ".zygosity" : zygosity})
	print ("The number of individuals who are " + description + " is " + str(cursor))

print ("--------------RESULTS--------------")
with open('Alleles.csv') as csvfile:
	linereader = csv.reader(csvfile, delimiter=',')
	for row in linereader:
		if row[0] == "Name":
			print (row[1])
		else:
			print (row[0])
			getCountOfMissing(getFindingsCollection('127.0.0.1', 27017), row[0])
			print ("*")
			getCountOfAlleles(getFindingsCollection('127.0.0.1', 27017), row[0], 0, row[1])
			getCountOfAlleles(getFindingsCollection('127.0.0.1', 27017), row[0], 1, row[1]+"/"+row[2])
			getCountOfAlleles(getFindingsCollection('127.0.0.1', 27017), row[0], 2, row[2])
			print ("\n")
