#!/bin/bash

# The following script is design for use at the NIH High Perfomance Computing resource Biowulf. Please adjust as necessary.

module load bbtools

#trim for quality and length
bbtools bbduk -Xmx24g in=/lscratch/$SLURM_JOBID/Sample.fastq out=/lscratch/$SLURM_JOBID/TRIMMED_Sample.fastq qtrim=rl trimq=20 minlen=25

#get rid of any singletons
bbtools repair in=/lscratch/$SLURM_JOBID/TRIMMED_Sample.fastq out=/lscratch/$SLURM_JOBID/Repaired_Trimmed_Sample.fastq outs=/lscratch/$SLURM_JOBID/Singletons_Trimmed_Sample.fastq repair

#map, ensure reference genome is in the working directory and that it has been previously indexed by bbmap.
bbtools bbmap in=/lscratch/$SLURM_JOBID/Repaired_Trimmed_Sample.fastq trd sam=1.3 ambiguous=toss out=/lscratch/$SLURM_JOBID/Sample.sam statsfile=/lscratch/$SLURM_JOBID/StatsBbmap_Sample.txt -Xmx28g

#load samtools package
module load samtools

#convert sam to bam
samtools view -bS /lscratch/$SLURM_JOBID/Sample.sam > /lscratch/$SLURM_JOBID/Sample.bam

# sort
samtools sort -o /lscratch/$SLURM_JOBID/Sample.sorted.bam -T /lscratch/$SLURM_JOBID/Sample.tmp /lscratch/$SLURM_JOBID/Sample.bam

#remove duplicates
samtools rmdup /lscratch/$SLURM_JOBID/Sample.sorted.bam /lscratch/$SLURM_JOBID/Sample.dedup.bam

#index
samtools index /lscratch/$SLURM_JOBID/Sample.dedup.bam

#rename the .dedup.bam.bai file, as required for Freebayes
cp /lscratch/$SLURM_JOBID/Sample.dedup.bam.bai /lscratch/$SLURM_JOBID/Sample.dedup.bai

# copy files that rylan and IgV will need back to your local storage
cp /lscratch/$SLURM_JOBID/Sample.dedup.bam.bai /data/user/DIRECTORY/.

cp /lscratch/$SLURM_JOBID/Sample.dedup.bai /data/user/DIRECTORY/.

cp /lscratch/$SLURM_JOBID/Sample.dedup.bam /data/user/DIRECTORY/.
