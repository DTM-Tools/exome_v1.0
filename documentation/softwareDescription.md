# DTM-Tools exome v1.0:Software Description

## Introduction
DTM-Tools was designed for prediction of blood antigen phenotypes from large NGS datasets. It is written in Python and designed to work in a Linux/Unix environment. It can only be run from the command line or through a Docker container – graphical and web interfaces will be implemented in the future.

The description and tutorial material presented in this document is intended as a complement, not a repetition, for the software’s publication. Please refer to ```Citation pending``` before reviewing this material.

## Source Code Structure

[Github main]: https://github.com/DTM-Tools/exome_v1.0

DTM-Tools is written in object-oriented notation as a series of stages (static methods) that create objects, execute, and return the output to the pipeline state. Logging levels are incorporated for future web server deployment.

The source code can be accessed through [Github][Github main]: https://github.com/DTM-Tools/exome_v1.0. Please refer to the DTM-Tools publication for additional details and laboratory validation data.

The following four subdirectories are included in exome v1.0:

_/build_: contains Docker definition files and the Docker deployment entrypoint.

_/databases_: contains three .csv files that define genomic coordinates and their interpretation rules. Described in further detail below.

_/rylantool_: contains DTM-Tools app python scripts

_/queryTools_: contains sample pymongo query scripts

## Dependencies

[pyvcf site]: https://pyvcf.readthedocs.io/en/latest/
[freebayes]:https://github.com/ekg/freebayes
[samtools]:http://www.htslib.org/

If run locally, DTM-Tools requires prior installation of [Samtools][samtools], [Freebayes][freebayes], and [pyvcf][pyvcf site]. 

However, a Docker container is available, which obviates the need for these installations and enhances reproducibility of DTM-Tools. Docker command lines can be found in _softwareTutorials.md_ in this /documentation directory.

You can read more about software reproducibility and the use of containers here:

[Peng paper]:https://www.ncbi.nlm.nih.gov/pubmed/?term=PMID%3A+22144613.

[Cito paper]:https://ieeexplore.ieee.org/document/7883438

[Piccolo paper]:https://www.ncbi.nlm.nih.gov/pubmed/?term=PMID%3A+27401684

[Schulz paper]:https://www.ncbi.nlm.nih.gov/pubmed/?term=PMID%3A+28163975

[Lobe paper]:https://www.ncbi.nlm.nih.gov/pubmed/?term=PMID%3A+27577463

[Wagholikar paper]:https://www.ncbi.nlm.nih.gov/pubmed/?term=PMID%3A+30012140

* Peng RD: [Reproducible research in computational science][Peng paper]. Science 2011; 334: 1226-7.
* Cito J, Gall HC: [Using Docker Containers to Improve Reproducibililty in Software Engineering Research][Cito paper];  2016 IEEE/ACM 38th IEEE International Conference on Software Engineering Companion. Austin, TX, 2016.
* Piccolo SR, Frampton MB: [Tools and techniques for computational reproducibility][Piccolo paper]. Gigascience 2016; 5: 30.
* Schulz WL, Durant TJ, Siddon AJ, et al.: [Use of application containers and workflows for genomic data analysis][Schulz paper]. J Pathol Inform 2016; 7: 53.
* Lobe M, Ganslandt T, Lotzmann L, et al.: [Simplified Deployment of Health Informatics Applications by Providing Docker Images][Lobe paper]. Stud Health Technol Inform 2016; 228: 643-7.
* Wagholikar KB, Dessai P, Sanz J, et al.: [Implementation of informatics for integrating biology and the bedside (i2b2) platform as Docker containers][Wagholikar paper]. BMC Med Inform Decis Mak 2018; 18: 66.

## Input File Format

* DTM-Tools takes .bam files and the associated index (.bai) file as input; these files should be stored in the directory path specified by the **‘-i’** option. Note that DTM-Tools currently employs Freebayes as its variant caller, which expects the index to have the ‘.bai’ extension, NOT ‘.bam.bai’. This requires renaming of the typical extension given by the samtools index command. Refer to our pre-processing script in this directory for sample commands for alignment and indexing. 

* DTM-Tools also requires the reference genome and its index as input files, stored in the directory path specified by the **‘-r’** command line option. Ensure that this is the exact same reference genome file employed when aligning the raw NGS FASTQ files. Be aware that the hg38 reference genome contains alternative haplotypes and patches which can lead to ambiguous alignments with some software. Note that in our sample preprocessing script we employ the bbmap aligner with the ‘ambiguous=toss’ option, which WILL result in loss of any reads derived from regions with alternative haplotypes (such as the _KEL_ gene). If using this script, limit the reference genome to the chromosomal assemblies only (no patches or haplotypes).

The reference genome version is hg19 by default. Hg38 can be specified with the **‘-g grch38’** option in the DTM-Tools command line. Currently, the expected names of the reference genome fasta files are hard-coded, you can modify them in rylantool/utils.py, getFastaFileName object (lines 19-26, shown below):

	def getFastaFileName(reference):
   	
   		if reference.lower() == "hg19":
        	return 'human_g1k_v37.fasta'
    	elif reference.lower() == "grch38":
        	return 'GRCh38.fasta'
    	else:
        	raise RylanParamException('Invalid reference type: %s' % reference)

DTM-Tools also has the option (**-a**) to report any additional genomic variants in genes of interest. For this purpose, the user is responsible for providing ANNOVAR output files. ANNOVAR output file extensions of ‘Annovar.variant\_function’ and ‘Annovar.exonic\_variant\_function’ are currently expected by the code. If the ANNOVAR output file name is different, you can modify the expected extension in rylantool/additional\_findings\_stage.py, compute\_annovar\_file\_names object (lines 31 and 32, shown below):


	def compute_annovar_file_names(self):
        base_name = self.rylan_params.get_patient_file()[:-10]
        input_dir = self.rylan_params.get_input_directory()
        file2_name = '%s/%s_Annovar.exonic_variant_function'%(input_dir, base_name)
        file1_name = '%s/%s_Annovar.variant_function'%(input_dir, base_name)
        return (file1_name, file2_name)

## DTM-Tools Database

The DTM-Tools command line also expects the genomic rules and interpretation databases stored in the directory path specified by the ‘**-r**’ command line option. 

DTM-Tools exome v1.0 database can be downloaded from the Github  _/database directory. Three files are expected: _ChromoList.csv, ChromoInDelList.csv, and Multi.csv_, each described individually below. 

[JSON dbsnpsite]: https://ncbiinsights.ncbi.nlm.nih.gov/2018/06/15/dbsnp-updates-json-refsnp-report-api

**Note that although downloaded files are in _.csv_ format, DTM-Tools reads them into memory as a JSON dictionary** . This will allow future versions of DTM-Tools to incorporate complex non-relational definitions, and renders the DTM-Tools backbone compatible with current [dbSNP JSON data]  [JSON dbsnpsite] file structure. 


* The user can open the donwloaded _.csv_ files and edit them in Excel, but must ensure that the line break remains as LF and not CRLF. 
* The first line of each _.csv_ database file is expected to be column headers and is ignored by the software. Note that entering allele information in the first line will not yield an error, however the information will not be interpreted in the output file.  
* DTM-Tools exome\_v1.0 database is based on the ISBT blood group alleles tables, but it does not use special characters to ensure compatibility with the DTM-Tools Python library. We provide detailed database nomenclature rules for the user in /documentation/databaseNomenclaturev1.pdf.
* Reference nucleotides must correspond exactly to the input reference genome. For example, although _SLC14A1 c.838A_ (_JK\*02_ or _JK\*B_) is listed as the reference in the ISBT tables, hg19 and hg38 lists _SLC14A1 c.838G_ (_JK\*01_ or _JK\*A_) as the reference at that position, and it is reflected as such in the DTM-Tools database.
* Note that nucleotide reference/variant classification and phenotype interpretation **must be provided in the PLUS strand**. For genes that are encoded in the negative strand, use the reverse complement.

### ChromoList.csv

ChromoList.csv is used to describe SNV interpretation rules and contains the following fields:

	A) Allele: unique code for each entry. See nomenclature rules for suggested nomenclature.
	
	B) Description: allele name per ISBT, without special characters as indicated in nomenclature rules. 
	
	C) Enabled: TRUE or FALSE
	
	D) System: blood group system number
	
	E) Gene: Name of gene
	
	F) Chromosome: chromosome number
	
	G) HG19Pos: nucleotide position of the SNV in hg19 (1-based)
	
	H) GRC38Pos: nucleotide position of the SNV in hg38 (1-based)
	
	I-L) A,C,G and T: enter ‘ref’ if it is considered the reference nucleotide (in the plus strand), enter ‘var’ if it is one of the documented variants, and enter ‘abn’ if it’s never been reported and you’d like DTM-Tools to flag it if found. Note that ‘ref’ is required, and only one ‘ref’ is allowed. More than one ‘var’ entry is possible.
	
	M-P) Aclass, Cclass,Gclass,Tclass: enter the phenotype classification you want DTM-Tools to assign to each of the nucleotides, if found in the plus strand. Remember to avoid any special characters that may interfere with the python code, such as asterisks or commas. Enter ‘-‘ if it is an ‘abn’ nucleotide call (i.e. unexpected nucleotide call for that genomic position). 
	
	Q) DP_REF: enter ‘x > #’ to denote the depth filter in this position if it is called homozygous reference
	
	R) MAPQ_REF: enter the minimum MAPQ filter in this position
	
	S) QUAL: qual value filter, applies only if there’s a variant nucleotide call
	
	T) DP: accepts minimum and maximum expected values for total depth, regardless of nucleotide call
	
	U) AO: filter for minimum number of detected alternate allele calls, applies only if there’s a variant nucleotide call
	
	V) AO/DP: filter for minimum fraction of alternate allele, applies only if there’s a variant nucleotide call
	
	W) Cell: ‘RBC’, ‘platelet’, or ‘PMN’
	
	X) Type: ‘Antigen’ vs ‘Enzyme’
	
	Y) Gensubtype: define if the SNV encodes for an epitope, a weakening variant, or a null change. Multiple values accepted in pipe-delimited format.
	
	Z-AC) Asubtype, Csubtype, Gsubtype, Tsubtype: specify what each nucleotide’s specific ‘Gensubtype’ (column Y) effect is.

### ChromoInDelList.csv

ChromoInDelList.csv is used to describe indel interpretation rules and contains the following fields:

	A) Allele: unique code for each entry. See nomenclature rules for suggested nomenclature.
	
	B) Description: allele name per ISBT, without special characters as indicated in nomenclature rules. 
	
	C) Enabled: TRUE or FALSE
	
	D) System: blood group system number
	
	E) Gene: Name of gene
	
	F) Chromosome: chromosome number
	
	G) HG19Start: Start nucleotide position in hg19 (see indel mapping note below)
	
	H) HG19End: End nucleotide position in hg19 (see indel mapping note below)
	
	I) GRC38Start: Start nucleotide position in hg38 (see indel mapping note below)
	
	J) GRC38End: End nucleotide position in hg38 (see indel mapping note below)
	
	K) Ref: reference nucleotide string expected in the plus strand
	
	L) Alt: alternate nucleotide string expected in the plus strand
	
	M) RefClass: phenotype classification for the ‘Ref’ nucleotide string (column K)
	
	N) AltClass: phenotype classification for the ‘Alt’ nucleotide string (column N)
	
	O) Type: ‘ins’ for insertion and ‘del’ for deletion; complex variations are not yet supported by DTM-Tools
	
	P) DP_REF: enter ‘x > #’ to denote the depth filter in this position if it is called homozygous reference
	
	Q) MAPQ_REF: enter the minimum MAPQ filter in this position for reference nucleotide strings
	
	R) MQM: filter for MQM value (mean mapping quality of observed alternate alleles)        
	
	S) MQMR: filter for MQMR value (mean mapping quality of observed reference alleles)
	
	T) QUAL: qual value filter, applies only if there’s a variant nucleotide call
	
	U) DP: accepts minimum and maximum expected values for total depth, regardless of nucleotide call
	
	V) AO: filter for minimum number of detected alternate allele calls, applies only if there’s a variant nucleotide call
	
	W) AO/DP: filter for minimum fraction of alternate allele, applies only if there’s a variant nucleotide call
	
	X) QR: filter for the sum of quality of the reference observations
	
	Y) QA: filter for the sum of quality of the alternate observations
	
	Z) Cell: ‘RBC’, ‘platelet’, or ‘PMN’
	
	AA) Function: ‘Antigen’ vs ‘Enzyme’
	
	AB) Gensubtype: define if the indel variant encodes for an epitope, a weakening variant, or a null change. 
	
	AC) Refsubtype: specific Gensubytpe for the reference nucleotide string.
	
	AD) Altsubtype: specific Gensubtype for the alternate nucleotide string.

#### Special considerations for Indels

DTM-Tools currently employs Freebayes, a haplotyper, for variant calling. Nucleotide coordinates for indels are 0-based start,1-based end. As a reference, the nucleotide coordinates displayed in the UCSC Genome Browser (graphical) are 1-based.

**Simulated insertion example:**

> For ATCGATCGATC  to ATCGATACGATC

Nucleotide A is added after the 6th nucleotide, T.

In DTM-Tools, if an insertion occurs after the Xth nucleotide, the start nucleotide coordinate is X-1, while the end nucleotide coordinate is X, so for this example:

> Start nucleotide coordinate specified in DTM-Tools database = 5
> 
> End nucleotide coordinate specified in DTM-Tools database = 6

**Simulated deletion example:**

> ATCGATACGATC to ATCGATCGATC

The nucleotide A between the 6th nucleotide, T, and the 8th nucleotide, C, is deleted.

In DTM-Tools, if a deletion occurs between the Xth and Yth nucleotide, the start nucleotide coordinate is X-1, while the end nucleotide coordinate is Y, so for this example:

> Start nucleotide coordinate in DTM-Tools database = 5 

> End nucleotide coordinate specified in DTM-Tools database = 8

### Multi.csv

Multi.csv is used for alleles defined by more than one SNV and/or indel, and contains the following fields:

	A) Allele: unique code for each entry. See nomenclature rules for suggested nomenclature.
	
	B) Description: allele name per ISBT, without special characters as indicated in nomenclature rules.
	
	C) List SNV: list of the ‘Allele’ unique codes (i.e. ‘Allele’ code for ChromoList.csv and ChromoInDel.csv) or each of the variants that define the allele
	
	D) RefClass: Phenotype classification if all the listed alleles in column C were not present in variant nucleotide form on the same haplotype
	
	E) AltClass: Phenotype classification if all the listed alleles in column C were found in alternate nucleotide form on the same haplotype

## DTM-Tools Output File
The DTM-Tools output file is in JSON format. It will be saved in the directory specified by the ‘**-o**’ option of the command line, and will have the name of the input file with the ‘_.DTM-Toolsout.json_’ extension. 

JSON format was selected, rather than a predicted individual phenotype output, to allow the user to perform powerful and complex queries in large cohorts using a non-relational database. This format will also allow for complete structural flexibility in the output of future software versions that incorporate copy number variations and genes with complex structural rearrangements. This format is also compatible with current genomic variation data file structures, such as [dbSNP] [JSON dbsnpsite].

* Please refer to the DTM-Tools publication for a description and illustration of the output file structure. 

* A sample output file is available here for download: /documentation/Sample.rylanout.json

* Sample non-relational database queries are available in /documentation/tutorials.md

* /documentation/tutorials.md also contains further tips on interpreting the DTM-Tools output JSON file

## Important Notes and Exceptions

* After a DTM-Tools run is completed, we recommend uploading JSON files to a non-relational database for querying. We recommend individual review of all cases whose classification is marked as “abnormality” and of all quality filter failures. 

* **Note that left-alignment and harmonization is not yet enabled in DTM-Tools**, which may result in variable representation of indels and ‘abnormality’ flags.

* Alleles defined by more than one SNV and/or indel (_Multi.csv_) are given a zygosity value and a certainty value in the JSON output file. Physical phasing is not enabled for DTM-Tools_exome, and thus phasing of heterozygous calls cannot be performed. 

Zygosity values are equivalent to those for SNVs and Indels: 
> 0 = homozygous reference

> 1 = heterozygous

> 2 = homozygous variant

Certainty values are:
 > 0 = uncertain (for example, three variants in heterozygous state)
 
 > 1 = certain (for example, two heterozygous variants and one homozygous result in a zygosity value of 1 and a certainty value of 1).
 
* DTM-Tools employs [Freebayes] [freebayes] for variant calling, which is haplotype-based. If an indel or a complex genomic variant is identified in a position listed in the _ChromoList.csv_ database (i.e. an expected SNV), the alternate nucleotide call will not match, and the allele will NOT be classified in the output database. This can be uncovered by running the DTM-Tools software with the ‘**-n**’ option, and inspection of the individual _.vcf_ file, which will be stored in the _/tmp_ directory. Thus, in the setting of exome sequencing, a missing allele in the output JSON can have 4 explanations: 

> 1) Indel or complex variant in an SNV-expected position

> 2) Trimming at pre-processing stage due to low quality/short sequence

> 3) Misalignment to paralogous gene

> 4) True homozygous deletion

**DTM-Tools is for research use only and is in continuous development**. Please contact the DTM-Tools developer at <celina.montemayorgarcia@nih.gov> for questions and to report any problems.