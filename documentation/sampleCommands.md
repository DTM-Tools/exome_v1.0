# DTM-Tools sample commands

## Running the Docker container

Running DTM-Tools as a Docker container is 
**highly recommended**; this eliminates the need to install any dependencies or activate environments. Refer to the DTM-Tools publication for a suplementary video tutorial illustrating the use of this Docker cotainer.

The following command, from the root directory, will build the container. Be sure to make any necessary modifications to the DTM-Tools source code before building the container, such as: 

* Expected names of reference genome file, and expected names for ANNOVAR output files: refer to the 'Input File Format' section in _/documentation/softwareDescription.md_ for further details. 
* Set the filter tresholds in _ChromoList.csv_ and _ChromoInDelList.csv_ that are appropriate for your dataset.
* If your reference genome names chromosomes as 'chr1', 'chr2', 'chr3', etc (rather than '1', '2', '3', etc), edit the source code accordingly. Refer to 'Input File Format' section in in _/documentation/softwareDescription.md_ for further details. 

Docker container build command:

	docker build -f build/rylan/Dockerfile -t dtmtool:latest . 
	
The following command will then run the container with the hg38 genome coordinates and six threads:

	docker run -v /DIRECTORY/PATH:/data -e INPUT_DIR=/data/Directory -e OUTPUT_DIR=/data/Directory -e PATIENT_FILE=Sample.dedup.bam -e NUMCORES=6 -e REF_DIR=/data/Directory -e REF_BUILD=grch38 -it dtmtool:latest

[freebayes]:https://github.com/ekg/freebayes
	
‘INPUT\_DIR’ should contain the file name given in the ‘PATIENT\_FILE’ environment variable, along with its index, and ANNOVAR output files. Remember that DTM-Tools employs [Freebayes][freebayes], which expects the index extension of _'.bai'_ rather than _'.bam.bai'_.

If you do not have ANNOVAR output files, add the following environment variable:

	-e RUN_ANNOVAR='False'

‘REF\_DIR’ should contain: ChromoList.csv, ChromoInDelList.csv, Multi.csv, genomic reference file and its index.

Hg19 is considered the default assembly value; if using hg19 the REF\_BUILD environment variable need not be specified.

## Running DTM-Tools locally

Activate the local environment (pyvcf) first. Make any necessary edits to the DTM-Tools source code and database (examples given in the previous section).

The following command line, from the root directory, executes DTM-Tools with hg38 database coordinates and six threads:

	python dtmtool/main.py -a -m 6 -r /DIRECTORY/PATH/REFERENCE -o /DIRECTORY/PATH/OUTPUT -i /DIRECTORY/PATH/INPUT -p Sample.bam -g grch38
	
/DIRECTORY/PATH/REFERENCE should contain: ChromoList.csv, ChromoInDelList.csv, Multi.csv, GRCh38.fasta, GRCh38.fasta.fai

/DIRECTORY/PATH/INPUT should contain the file name given in the **-p** option, along with its index, and ANNOVAR output files if **-a** is employed.

To keep a copy of each generated .vcf file add a **-n** option to the command line.

If you don’t have ANNOVAR output files, omit the **-a** option.

Hg19 is considered the default assembly build value, in this case it does not need to be specified and the following command line will work:

	python dtmtool/main.py -a -m 6 -r /DIRECTORY/PATH/REFERENCE -o /DIRECTORY/PATH/OUTPUT -i /DIRECTORY/PATH/INPUT -p Sample.bam
	
## Evaluating the JSON output file

The JSON output file and its rationale is also discussed in /_documentation/softwareDescription.md_ - please review that document first. **Figure 2** the DTM-Tools publication provides an illustration and description of each field in the DTM-Tools JSON output.

To view the JSON output file in a visually-friendly format, use any online JSON formatter (search for "JSON formatter" in any web search engine). Some examples are:

[sample1]:https://jsonformatter.curiousconcept.com/

[sample2]:https://www.jsonformatter.io/

[sample3]:https://jsonformatter-online.com/

[Open formatter example 1][sample1]

[Open formatter example 2][sample2]

[Open formatter example 3][sample3]

The user can then download the formatted JSON file and visualize with any text editor. In addition, the Firefox web browser will automatically open any _.json_ file in an expandible/collapsible format.

A **supplementary video tutorial** ```hyperlink pending``` is available in th DTM-Tools publication, illustrating the structure and interpretation of the JSON output file, and the use of pymongo script to query individual and cohort results.

	
## Pymongo scripts

Pymongo scripts are provided in the /queryTools directory. These are provided as a guide; we upload our JSON files into MongoDB for querying and analysis but any non-relational database could be used. A ***video tutorial*** is provided in the DTM-Tools publication supplementary material and illustrates the use of these scripts ```hyperlink pending```.

[pymongo site]:https://api.mongodb.com/python/current/

DTM-Tools exome_v1.0 release contains the following scripts which require [pymongo][pymongo site]. Note that authentication is not encoded in these scripts; **authentication is the responsibility of the user.**

### AlleleCountMissingFiltered.py

This script provides an output that indicates, per row entry in each of the three databases, how many documents are classified as homozygous reference (i.e. zygosity code = 0), heterozygous (zygosity code =1) or homozygous variant (zygosity code =2). It also indicates how many documents were flagged as ‘filtered’ due to a failed quality filter, and how many documents did not contain such allele entry (missing in the input .bam file). Requires Python 3.x, and expects _Alelles.csv_ in the running directory. Note that some genomic positions have more than one variant nucleotide classification, in that case this script will not specify; use direct database queries to dissect these cases.

The collection and table names are specified in line 10, remember to edit this before running the script:

	return client.collectionName.tableName

Command line to run the script:

	python AlleleCountMissingFiletered.py > /DIRECTORY/PATH/outputFile.txt

### MultiAlleleCount.py

This script is a variation of _AlleleCountMissingFiltered.py_ that addresses alleles that are defined by more than one SNV/Indel. Requires _AllelesMultiKKD.csv_ in the running directory.

Collection and table names are specified in line 10, as above. 
	
### AddVariantsPerGene.py

This script provides a breakdown of the genomic variants identified, per ANNOVAR, whose coordinates are NOT described in the DTM-Tools database used in that particular run. 

Output provides a count of total exonic, intronic, 3’UTR, and 5’UTR variants per gene. Note that gene names are determined by ANNOVAR and may vary. For exonic variants, the script provides the breakdown count for nonsynonymous SNV, synonymous SNV, nonframeshifts (substitutions, insertions, and deletions), frameshifts (substitutions, insertions, and deletions), stop gains, stop losses, and unknowns - as defined by ANNOVAR. Requires Python 3.x. 

The collection and table names are specified in line 10, remember to edit this before running the script:

	return client.collectionName.tableName

Command line:

	python AddVariantsPerGene.py > /DIRECTORY/PATH/outputFile.txt
	
### AddVariantsInputGene.py
This script is a variation of _AddVariantsPerGene.py_. It allows the user to specify the gene names to query for novel variants, and it has the option to report intronic changes per ANNOVAR.

Collection and table names are specified in line 10, as above.

Note that the section to report intron variants, which is very time-consuming, is commented out in the tracked version of this script. To include it, simply remove the 3 tick marks ```'''```at the beginning and at the end of this section:

    '''idListIntronic=[]
    countIntronic=0
    for id in alleles.keys():
        if id not in idListIntronic:
            if alleles[id]['location']=="intronic":
                countIntronic = countIntronic + 1
                idListIntronic.append(id)
    print ('Total intronic changes:',+countIntronic)'''
    
 Deafult gene names to scan are KEL, SLC14A1, and ACKR1; these can be edited in line 21:
 
 	geneList=['KEL','SLC14A1','ACKR1']

Command line:

	python AddVariantsInputGene.py > /DIRECTORY/PATH/outputFile.txt

	
## Sample MongoDB queries

Any non-relational database can be employed. We provide here some commands that the development team uses to import and analyze documents in MongoDB for your reference.

To import in bulk all .rylanout.jsons (DTM-Tools output files) saved in a directory:

	find . -regex '.*/[^/]*.json' | xargs -L 1 mongoimport --db databaseName --collection tableName --file
	
To list the runID and Input File name for a document that contains a specific additionalFinding:

	db.databaseName.find({"additionalFindings.CHROM19:45321793:G:T.gene":"BCAM"},{"_id":0,"runId":1,"inputFile":1})

	
To count the documents that predict homozygosity for s8a (Fyb) and homozygosity for the s8b (GATA box mutation):

	db.databaseName.count({$and:[{"findings.s8a.zygosity":2},{"findings.s8b.zygosity":2}]})
	
To list the runIDs for all documents that predict homozygosity for s8a (Fyb) and homozygosity for the s8b (GATA box mutation):

	db.databaseName.find({$and:[{"findings.s8a.zygosity":2},{"findings.s8b.zygosity":2}]},{"_id" : 0,"runId":1})
	
To list the runIDs for all documents that have a certain allele (s6a in this case) not classified and not filtered:

	db.databaseName.find({$nor:[{"findings.s6a.determination":"classified"},{"findings.s6a.determination":"filtered"}]},{"_id" : 0,"runId":1})
	
**DTM-Tools is for research use only and is in continuous development**. Please contact the DTM-Tools developer at <celina.montemayorgarcia@nih.gov> for questions and to report any problems.