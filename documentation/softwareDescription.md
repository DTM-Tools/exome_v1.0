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

[Peng paper]:hhttps://www.ncbi.nlm.nih.gov/pubmed/?term=PMID%3A+22144613.

[Cito paper]: https://ieeexplore.ieee.org/document/7883438

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

* The user can open these files and edit them in Excel, but must ensure that the line break remains as LF and not CRLF. 
* The first line of each .csv database file is expected to be column headers and is ignored by the software. Note that entering allele information in the first line will not yield an error, however the information will not be interpreted in the output file.  
* DTM-Tools exome v1.0 database is based on the ISBT blood group alleles tables, but it does not use special characters to ensure compatibility with the DTM-Tools Python library. We provide detailed database nomenclature rules for the user in /documentation/databaseNomenclaturev1.pdf.
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
