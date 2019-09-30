from pymongo import MongoClient
import pymongo
import sys
import csv

def getFindingsCollection(ipAddr, srvport): #Function calls specific database and collection from Mongo, requires IP and server numbers as inputs
	client = MongoClient(host=ipAddr, port=srvport)
	db = client.rylanFindings #After client, write out database name that will be queried
	#return db.EvenMoreResults #Returns specific collection name and "pathway", anywhere database is needed, function can be used to call database
	return client.collectionName.tableName

def getCountOfMissing(findingscollection, allele):
	cursor = findingscollection.count({"findings." + allele + ".determination" : "MISSING"})
	print ("The number of individuals with this locus MISSING in input file is " + str(cursor))
	cursor = findingscollection.count({"findings." + allele + ".determination" : "filtered"})
	print ("The number of individuals with this locus FILTERED is " + str(cursor))

def getCountOfAlleles(findingscollection, allele, zygosity, description): #Function calls specific number of matched alleles, prints number as well as specific patient IDs for unique variants
	cursor = findingscollection.count({"findings." + allele + ".zygosity" : zygosity}) #cursor contians all number of patients with pathway "findings.AlleleName.zygosity.zygosityValue"
	print ("The number of individuals who are " + description + " is " + str(cursor))

print ("--------------RESULTS--------------")
with open('Alleles.csv') as csvfile: #Accesses the excel file as titled, states that row[0] is allele name, row[1] is zygosity value, row[2] is description, and rows[3] and [4] describe whether the patientID will be printed
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
