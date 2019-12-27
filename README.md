# DTM-Tools classifier
This is the DTM-Tools source code page. 

DTMTools was created by the Department of Transfusion Medicine of the National Institutes of Health, in loving memory of our patient Rylan Jaxson Macomb.

You can learn more about DTMTools in the _/documentation_ directory. **We strongly encourage users to review all documentation material before downloading and running the DTM-Tools code**.

## Download Reference Genome

DTM-Tools currently supports hg19 and hg38. Note that hg38 reference downloads often contain alternative haplotypes and patches, which may lead to ambiguous mapping. 

Currently, DTM-Tools expects the following exact filenames for hg19 reference build and index:

* human\_g1k\_v37.fasta
* human\_g1k\_v37.fasta.fai

For hg38, expected file names are:

* GRCh38.fasta
* GRCh38.fasta.fai

Expected reference file names can be modified in _rylantool/utils.py_.

## Dependencies
Please install the following tools before running:

* samtools
* freebayes

Ensure that your Linux system has the development packages installed, including support for **make**.

DTM-Tool's output is a series of JSON files that can be subsequently queried an anlalyzed with a non-relational database. We employ **MongoDB** and provide our analysis scripts in the _querytools_ directory.

## Preparing Aligned (.bam) Files

DTM-Tools takes sorted bam files as input, along with the corresponding index. Currently DTM-Tools employs Freebayes for variant calling; validation was performed with bbmap alignment. The  _/documentation_ directory provides a sample trimming and aligning pipeline. 

Example sorting and indexing pipeline for an aligned file:	

	# Sort the BAM file
	samtools sort -o FILENAME.mapped.sorted.bam -T /tmp/FILENAME.mapped.aligned.tmp FILENAME.mapped.bam
	# Remove duplicates and index output
	samtools rmdup HG00264.mapped.sorted.bam HG00264.mapped.dedup.bam
	samtools index HG00264.mapped.dedup.bam 
	# DTM-Tools expects .bam and .bai to have the same name (change extension)
	mv HG00264.mapped.dedup.bam.bai HG00264.mapped.dedup.bai
	
## Genomic Interpretation Rules Database

DTM-Tools interprets Freebayes vcf results in relation to an interpretation rules database. For blood group interpretation, we provide databases in .csv format (_ChromoList.csv, ChromoInDelList.csv and Multi.csv_) located in _/databases_. Please refer to the _/documentation_ directory and to our publication for further details.

## Additional Genomic Variants in Genes of Interest
[Annovar webpage]: http://annovar.openbioinformatics.org/en/latest/
DTM-Tools accepts [ANNOVAR][Annovar webpage] output files to annotate additional findings in the genes of interest (i.e. genomic variants not present in the DTM-Tools database). Users must run Annovar anlaysis independently and store the corresponding output files in the input directory.

## Running the Classifier

To run the DTMTools classifier, you will need to set up Python 3.x. We recommend that you use a virtual environment to configure your running setup, with the following dependencies:

	pip install pyvcf
	
To run DTM-Tools exome v1.0, this is the command line and current supported options:

	python rylantool/main.py [-a] [-n] -m <num_cores> -r <refdir> -o <outputdir> -i <inputdir> -p <patientfile.dedup.bam> -g <hg19|grch38>
		
Where:

- **-a**: Run annovar analysis to capture additional variants in the blood group genes. Annovar must be run independently by the user; DTM-Tools expects Annovar output files in the input directory. Expected file names are specified in _rylantools/additional\_findings\_stage.py_.
- **-n**: Do not clean up temporary files (for debugging).
- **-m**: Number of parallel processes to make the analysis faster. For use in multicore machines. Recommended one per core. Maximum number of threads = number of chromosomes in the database.
- **-r**: Directory where the FASTA genomic reference files and respective index are found. Software also expects genomic coordinate .csv files in this directory (_ChromoList.csv, ChromoInDelList.csv, Multi.csv_).
- **-o**: Ouptut will be written as a .json file in this specified directory.
- **-i**: Directory where aligned, sorted bam files and indexes are located. If _-a_ option is used, DTM-Tools expects Annovar output files in this directory.
- **-p**: Exact name of the .bam file to be anlayzed in this run. Must be located in _inputdir_ (specified by the **-i** option).
- **-g**: Genome build coordinates to employ in the interpretation;  *hg19* or *grch38* values are supported today (case sensitive). Default value is _hg19_.

To display the usage line:

```
python rylantool/main.py -h
```

## Building Docker Container

The following command will build a docker container, if run from the root of this repository:

```
docker build -f build/rylan/Dockerfile -t dtmtool:latest .
```

## Running DTM-Tools in Docker
To run the container in automatic mode:
```
docker run -v <host_dir>:<container_dir> -e <env1>='<val1>' ... -it dtmtool:latest
```

Environment variables are the only way to pass environment to the containers, via the '-e' switch (repeat once per variable).
The supported variables and default values are:

```
* RYLAN_INPUT_DIR
* RYLAN_OUTPUT_DIR
* RYLAN_REF_DIR
* RYLAN_PATIENT_FILE
* RYLAN_PATIENT_ID = 'unspecified'
* RYLAN_REF_BUILD = 'hg19'
* RYLAN_NUMCORES = '2'
* RUN_ANNOVAR = 'True'
```
## Legal

This software/database is “United States Government Work” under the terms of the United States Copyright Act. It was written as part of the authors’ official duties for the United States Government and thus cannot be copyrighted. The public may use the software/database on a worldwide and royalty-free basis for any purpose and may reproduce and prepare derivative works without limitation. Although all reasonable efforts have been taken to ensure the accuracy and reliability of the software/database and associated data, the National Institutes of Health (NIH), the NIH Clinical Center, and the U.S. Government disclaim all warranties as to performance, merchantability or fitness for any particular purpose.

This software is intended for research use only.

In any work or product derived from this material, proper attribution of the authors as the source of the software/database or data should be made, using the following citation:

```
*Citation pending
```


Questions? Contact the DTM-Tools developer at <celina.montemayorgarcia@nih.gov>. 
We welcome all comments and contributions.

