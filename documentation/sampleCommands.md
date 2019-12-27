# DTM-Tools sample commands

## Running the Docker container

Running RyLAN as a Docker container is 
**highly recommended**; this eliminates the need to install any dependencies or activate environments.

The following command, from the root directory, will build the container. Be sure to make any necessary modifications to _/rylantools_ before building the container (such as expected names of reference genome file, and expected names for ANNOVAR output files; refer to the 'Input File Format' section in _/documentation/softwareDescription.md_ for further details).

	docker build -f build/rylan/Dockerfile -t rylantool:latest 
	
The following command will then run the container with the hg38 genome coordinates and six threads:

	docker run -v /DIRECTORY/PATH:/data -e RYLAN_INPUT_DIR=/data/Directory -e RYLAN_OUTPUT_DIR=/data/Directory -e RYLAN_PATIENT_FILE=Sample.dedup.bam -e RYLAN_NUMCORES=6 -e RYLAN_REF_DIR=/data/Directory -e RYLAN_REF_BUILD=grch38 -it rylantool:latest
	
‘RYLAN\_INPUT\_DIR’ should contain the file name given in the ‘RYLAN\_PATIENT\_FILE’ environment variable, along with its index, and ANNOVAR output files.

If you do not have ANNOVAR output files, add the following environment variable:

	-e RUN_ANNOVAR='False'

‘RYLAN\_REF\_DIR’ should contain: ChromoList.csv, ChromoInDelList.csv, Multi.csv, genomic reference file and its index.

* Reference genome names are currently hard coded. Refer to the 'Input File Format' section in /documentation/softwareDescription.md for instructions on editing the expected file name.
* Hg19 is considered the default assembly value.

## Running DTM-Tools locally

Activate the local environment (pyvcf) first. The following command line, from the root directory, executes DTM-Tools with hg38 database coordinates and six threads:

	python rylantool/main.py -a -m 6 -r /DIRECTORY/PATH/REFERENCE -o /DIRECTORY/PATH/OUTPUT -i /DIRECTORY/PATH/INPUT -p Sample.bam -g grch38
	
/DIRECTORY/PATH/REFERENCE should contain: ChromoList.csv, ChromoInDelList.csv, Multi.csv, GRCh38.fasta, GRCh38.fasta.fai

* Reference genome names are currently hard coded. Refer to the 'Input File Format' section in /documentation/softwareDescription.md for instructions on editing the expected file name.

•	/DIRECTORY/PATH/INPUT should contain the file name given in the -p option, along with its index, and ANNOVAR output files if -a is employed.

To keep a copy of each generated .vcf file add a -n option to the command line.

If you don’t have ANNOVAR output files, omit the -a option.

Hg19 is considered the default assembly build value, in this case the following command line will work:

	python rylantool/main.py -a -m 6 -r /DIRECTORY/PATH/REFERENCE -o /DIRECTORY/PATH/OUTPUT -i /DIRECTORY/PATH/INPUT -p Sample.bam
	
## Understanding the JSON output file

The JSON output file and its rationale is discussed in /documentation/softwareDescription.md - please review first.

To view the .JSON output file in a visually-friendly format, use any online JSON formatter (look for "JSON formatter in any web search engine). Some examples are:

[formatter1]: [https://jsonformatter.curiousconcept.com/]

[formatter2]: [https://www.jsonformatter.io/]

[formatter3]: [https://jsonformatter-online.com/]

[Open formatter example 1][formatter1]

[Open formatter example 2][formatter2]

[Open formatter example 3][formatter3]

The user can then download the formatted JSON file and visualize with any text editor. In addition, the Firefox web browser will automatically open any _.json_ file in an expandible/collapsible format.

* Understanding the structure of the JSON file is paramount to be able to perform the proper final phenotype prediction, and to design complex queries for large cohorts. Please refer to the DTM-Tools publication for an illustration and description of each field in the DTM-Tools JSON output. A sample output file is available here for download: /documentation/Sample.rylanout.json.

	
## Pymongo scripts

Pymongo scripts are provided in the /queryTools directory. These are provided as a guide, we upload our JSON files into MongoDB for querying and analysis but any non-relational database could be used.

[pymongo site]:https://api.mongodb.com/python/current/

DTM-Tools exome_v1.0 release contains the following two scripts which require [pymongo][pymongo site]:

### AlleleCountMissingFiltered.py

This script provides an output that indicates, per row entry in each of the three databases, how many documents are classified as homozygous reference (ie zygosity code = 0), heterozygous (zygosity code =1) or homozygous variant (zygosity code =2). It also indicates how many documents were flagged as ‘filtered’ due to a failed quality filter, and how many documents did not contain such allele entry (missing in the input .bam file). Requires Python 3.x. Note that authentication is not encoded in this instance; authentication is the responsibility of the user.

The collection and table names are specified in line 10, remember to edit this before running the script:

	return client.collectionName.tableName

Command line to run the script:

	python AlleleCountMissingFiletered.py > /DIRECTORY/PATH/outputFile.txt
	
### AddVariantsPerGene.py

This script provides a breakdown of the genomic variants identified, per ANNOVAR, whose coordinates are NOT described in the DTM-Tools database used in that particular run. 

Output provides a count of total exonic, intronic, 3’UTR, and 5’UTR variants per gene. Note that gene names are determined by ANNOVAR and may vary. For exonic variants, the script provides the breakdown count for nonsynonymous SNV, synonymous SNV, nonframeshifts (substitutions, insertions, and deletions), frameshifts (substitutions, insertions, and deletions), stop gains, stop losses, and unknowns - as defined by ANNOVAR. Requires Python 3.x. Note that authentication is not encoded in this instance; authentication is the responsibility of the user.

The collection and table names are specified in line 10, remember to edit this before running the script:

	return client.collectionName.tableName

Command line:

	python AddVariantsPerGene.py > /DIRECTORY/PATH/outputFile.txt
	
## Sample MongoDB queries

Any non-relational database can be employed. We provide here some commands that the development team uses to import, export and analyze documents in MongoDB for your reference.

To import in bulk all .rylanout.jsons in a directory:

	find . -regex '.*/[^/]*.json' | xargs -L 1 mongoimport --db databaseName --collection tableName --file
	
To list the runID and Input File name for a document that contains a specific additionalFinding:

	db.databaseName.find({"additionalFindings.CHROM19:45321793:G:T.gene":"BCAM"},{"_id":0,"runId":1,"inputFile":1})

	
To count the documents that predict homozygosity for sa8 (Fyb) and homozygosity for the s8b (GATA box mutation):

	db.databaseName.count({$and:[{"findings.sa8.zygosity":2},{"findings.sb8.zygosity":2}]})
	
To list the runIDs for all documents that predict homozygosity for sa8 (Fyb) and homozygosity for the s8b (GATA box mutation):

	db.databaseName.find({$and:[{"findings.sa8.zygosity":2},{"findings.sb8.zygosity":2}]},{"_id" : 0,"runId":1})
	
To list the runIDs for all documents that have a certain allele (s6a in this case) not classified and not filtered:

	db.databaseName.find({$nor:[{"findings.s6a.determination":"classified"},{"findings.s6a.determination":"filtered"}]},{"_id" : 0,"runId":1})