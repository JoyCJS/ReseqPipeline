#!/bin/bash

#Author: Joy Sung

#Codes below were used specifically for the SNP analyses of all 48 samples in this study.
#Before running this script, make a new directory to save the required files listed below:
#96 fastq.gz files
#Reference.fna
#Names.txt

#Unzip the fastq.gz files.
gunzip *.gz

#Make a list of paired sample names so samples can be run using a loop.
ls *.fastq | wc -l #This should print the number of fastq files.
ls *.fastq > a.txt
#Based on the fastq file names, modify the list to shorter names by removing the tails.
#Example file names: "LBK22-A01_S1_L001_R1_001.fastq", "LBK22-A01_S1_L001_R2_001.fastq", "LBK22-A06_S41_L001_R1_001.fastq", and "LBK22-A06_S41_L001_R2_001.fastq".
#Example shortened paired sample names: "LBK22-A01_S1" and "LBK22-A06_S41".
sed -i 's/_L001_R1_001.fastq//g' a.txt #Modify this part based on the fastq file names.
sed -i 's/_L001_R2_001.fastq//g' a.txt #Same here.
sort a.txt | uniq > names.txt
wc names.txt #This should print the number of samples (half of the number of fastq files).
rm a.txt

#Create a new directory for analysis.
mkdir GATK
cd GATK/
#Move required files into this directory.
mv ../names.txt .
mv ../Names.txt .
mv ../Reference.fna .

#The reference used in this study is an artificial reference of 7,755 sequences (used as 7,755 chromosomes; including 2,770 targets and 4,985 corresponding paralogs/homoeologs).
#Reference bwa, create file index, and generate sequence dictionary.
bwa index -a bwtsw Reference.fna
samtools faidx Reference.fna
mv Reference.fna Reference.fasta
java -jar /opt/picardtools/picard-tools-1.98/CreateSequenceDictionary.jar REFERENCE=Reference.fasta OUTPUT=Reference.dict
mv Reference.fasta Reference.fna
#There should be 8 files (Reference.dict, Reference.fna, Reference.fna.amb, Reference.fna.ann, Reference.fna.bwt, Reference.fna.fai, Reference.fna.pac, and Reference.fna.sa).

#Run bwa; modify this part based on your file names if needed.
while read A; do bwa aln -t 12 Reference.fna ../"$A"_L001_R1_001.fastq > $A.1.bwa; done < names.txt
while read A; do bwa aln -t 12 Reference.fna ../"$A"_L001_R2_001.fastq > $A.2.bwa; done < names.txt
ls *.bwa | wc -l #This should print the number of fastq files.

#Run sampe; modify this part based on your file names if needed.
while read A; do bwa sampe Reference.fna $A.1.bwa $A.2.bwa ../"$A"_L001_R1_001.fastq ../"$A"_L001_R2_001.fastq > $A.sam; done < names.txt
ls *.sam | wc -l #This should print the number of samples.

#Run bam files.
while read A; do samtools view -F 4 -Sbh $A.sam > $A.bam; done < names.txt
ls *.bam | wc -l #This should print the number of samples.

#Run picardtools with RGPU information.
while read A B C D; do java -jar /opt/picardtools/picard-tools-1.98/AddOrReplaceReadGroups.jar I=$A.bam O=$A.RG.bam RGID=$B RGLB=$C RGPL=illumina RGPU=$D RGSM=$C SORT_ORDER=coordinate CREATE_INDEX=true; done < Names.txt
ls *.bam | wc -l #This should print the number of fastq files.
ls *.bai | wc -l #This should print the number of samples.
awk '{print $1}' Names.txt > NamesPE.txt

#SortSam.
while read A; do java -jar /opt/picardtools/picard-tools-1.98/SortSam.jar I=$A.RG.bam O=$A.RG.sorted.bam SORT_ORDER=coordinate; done < NamesPE.txt
while read A; do java -jar /opt/picardtools/picard-tools-1.98/MarkDuplicates.jar I=$A.RG.sorted.bam O=$A.RG.sorted.dedup.bam M=$A.RG.sorted.bam.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000; done < NamesPE.txt
while read A; do java -jar /opt/picardtools/picard-tools-1.98/BuildBamIndex.jar I=$A.RG.sorted.dedup.bam; done < NamesPE.txt

#Run GATK to identify SNPs against reference for each sample using Haplotype Caller with parameter for 0/0 genotype (reference confidence).
#For more information, please check the GATK Guidebook (Page 379/400).
#Also check this paper "Computational Exome and Genome Analysis" by Peter N. Robinson, Rosario Michael Piro, and Marten Jager (2017) published at CRC Press (Chapter 13.10, Page 199).
while read A; do java -jar /opt/GenomeAnalysisTK/GenomeAnalysisTK.jar -R Reference.fna -T HaplotypeCaller -variant_index_type LINEAR -variant_index_parameter 128000 --emitRefConfidence GVCF -I $A.RG.sorted.dedup.bam -stand_call_conf 20 -stand_emit_conf 20 --allow_potentially_misencoded_quality_scores -o $A.RG.raw.snps.indels.vcf; done < NamesPE.txt
#Time varies but it took almost one day for us to finish running 48 samples.
ls *.vcf | wc -l #This should print the number of samples.
ls *.vcf.idx | wc -l #This should print the number of samples.

#Identify SNPs among samples.
while read A; do /opt/tabix/bgzip $A.RG.raw.snps.indels.vcf; done < NamesPE.txt
while read A; do /opt/tabix/tabix $A.RG.raw.snps.indels.vcf.gz; done < NamesPE.txt
ls *.vcf.gz | wc -l #This should print the number of samples.
ls *.vcf.gz.tbi | wc -l #This should print the number of samples.
#Vcf merge.
export PERL5LIB=/opt/vcftools/vcftools_0.1.11/perl
export PATH=$PATH:/opt/tabix/
ls *.RG.raw.snps.indels.vcf.gz | wc -l #This should print the number of samples.
vcf-merge *.RG.raw.snps.indels.vcf.gz | /opt/tabix/bgzip -c > pooledHC.vcf.gz
#Time varies but it took around 8 hours for us to finish running 48 samples.

#Unzip the pooled vcf.gz file.
gunzip pooledHC.vcf.gz
#The output file "pooledHC.vcf" can be downloaded for further downstream analyses.
#Reminder: Not all of the SNPs are targets; in this study, the targeted SNPs are at the 151st positions on the target sequences only.
#The next step is to filter the vcf file; you can use any other filtering method, or follow the codes that we used in this study (see "2_Data_Filtering.sh").