from pymongo import MongoClient
import pymongo
import sys
import csv
import statistics


def getAlleleList(File1,File2):
    idList=[]
    with open(File1) as csvfile:
        linereader = csv.reader(csvfile,delimiter=',')
        for row in linereader:
            if (row[0]) != 'Allele':
                idList.append(row[0])
    with open(File2) as csvfile:
        linereader2 = csv.reader(csvfile,delimiter=',')
        for row2 in linereader2:
            if (row2[0]) != 'Allele':
                idList.append(row2[0])
    return idList

def getFindingsCollection(ipAddr, srvport):
	client = MongoClient(host=ipAddr, port=srvport)
	db = client.collectionName
	return client.collectionName.tableName

def getAllDocsTree(findingscollection, alleles, metric):

    docs = findingscollection.find()
    runs = dict()
    for document in docs:
        alleledet = dict()
        for allele in alleles:
            if (document['findings'].get(allele, None) is not None):
                captured_findings = dict()
                if (document['findings'][allele].get(metric, None) is not None):
                    captured_findings[metric] = document['findings'][allele][metric]
                alleledet[allele] = captured_findings
        runs[document['runId']] = alleledet

    return runs

def flattenRunsMap(all_runs):
    flat_all_runs = dict()
    for runId in all_runs.keys():
        for allele in all_runs[runId].keys():
            flat_key = runId + "." + allele
            flat_all_runs[flat_key] = all_runs[runId][allele]
    return flat_all_runs

def main(argv):
    metric = argv[0]
    metric_name = argv[1]
    alleles = getAlleleList('ChromoList.csv','ChromoInDelList.csv')
    collection = getFindingsCollection('127.0.0.1', 27017)
    all_runs = getAllDocsTree(collection, alleles, metric)
    all_runs_flat = flattenRunsMap(all_runs)

    flat_keys = all_runs_flat.keys()


    print ('%s stats per allele:'%(metric_name))
    GenMAPQ=[]
    for allele in alleles:
        flat_keys_allele = list(filter(lambda x: x.endswith('.' + allele) == True, flat_keys))
        count_allele=0
        acc_allele=[]
        for a_key in flat_keys_allele:

            if (all_runs_flat[a_key].get(metric, None) is not None):
                acc_allele.append(all_runs_flat[a_key][metric])
                count_allele = count_allele + 1
                GenMAPQ.append(all_runs_flat[a_key][metric])
        print ('\n')
        print ('Stats for '+allele +':')
        if count_allele != 0:
            print ('The mean is: ' + str (statistics.mean(acc_allele)))
            print ('The minimum is: ' + str(min(acc_allele)))
            print ('The max is: ' + str(max(acc_allele)))
            if len(acc_allele) > 1:
                print ('The STDEV is: ' + str(statistics.stdev(acc_allele)))
        else:
            print ('No entries found')


    print('\n')
    print ('Overall %s stats:'%(metric_name))
    print ('The mean is: ' + str (statistics.mean(GenMAPQ)))
    print ('The minimum is: ' + str(min(GenMAPQ)))
    print ('The max is: ' + str(max(GenMAPQ)))
    if len(acc_allele) > 1:
        print ('The STDEV is: ' + str(statistics.stdev(GenMAPQ)))

if __name__ == "__main__":
   main(sys.argv[1:])
