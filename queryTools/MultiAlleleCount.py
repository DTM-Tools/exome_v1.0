from pymongo import MongoClient
import pymongo
import sys
import csv

def getFindingsCollection(ipAddr, srvport):
	client = MongoClient(host=ipAddr, port=srvport)
	db = client.collectionName
	return client.collectionName.tableName

def getCountOfMissing(findingscollection, allele):
	cursor = findingscollection.count({"findings." + allele + ".determination" : "not classfied"})
	print ("The number of individuals with this allele NOT CLASSIFIED is " + str(cursor))
	cursor = findingscollection.count({"findings." + allele + ".determination" : "filtered"})
	print ("The number of individuals with allele FILTERED is " + str(cursor))

def getCountOfAlleles(findingscollection, allele, zygosity, description):
	cursor = findingscollection.count({"findings." + allele + ".zygosity" : zygosity})
	print ("The number of individuals who are " + description + " is " + str(cursor))

print ("--------------RESULTS--------------")
with open('AllelesMultiKKD.csv') as csvfile:
	linereader = csv.reader(csvfile, delimiter=',')
	for row in linereader:
		if row[0] == "Name":
			print (row[1])
			print ("PLEASE REVIEW CERTAINTY VALUES")
			print("\n")
		else:
			print (row[0])
			print (row[3])
			getCountOfMissing(getFindingsCollection('127.0.0.1', 27017), row[0])
			print ("*")
			getCountOfAlleles(getFindingsCollection('127.0.0.1', 27017), row[0], 0, row[1])
			getCountOfAlleles(getFindingsCollection('127.0.0.1', 27017), row[0], 1, row[1]+"/"+row[2])
			getCountOfAlleles(getFindingsCollection('127.0.0.1', 27017), row[0], 2, row[2])
			print ("\n")
