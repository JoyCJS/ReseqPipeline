#!/bin/bash

#Author: Joy Sung

#Continued after the "1_Variant_Calling.sh" steps.
#Codes below were used specifically for the SNP analyses of all 48 samples in this study.
#To run this script, save the required file listed below:
#2770TargetInfo.txt

#Reminder: Not all of the SNPs are targets; in this study, the targeted SNPs are at the 151st positions on the target sequences only.
#The next step is to filter the vcf file; you can use any other filtering method, or follow the codes that we used in this study (see below).

#Filter by GQ and DP (typically GQ < 30 and DP < 10, or with loose threshold GQ < 6 and DP < 2).
#I followed the steps from the bottom of this webpage: https://evodify.com/gatk-in-non-model-organism/.
java -jar /opt/GenomeAnalysisTK/GenomeAnalysisTK.jar -T VariantFiltration -R Reference.fna -V pooledHC.vcf -G_filter "GQ < 30 || DP < 10" -G_filterName "Filtered" -o FilteredStep1.vcf
java -jar /opt/GenomeAnalysisTK/GenomeAnalysisTK.jar -T SelectVariants -R Reference.fna -V FilteredStep1.vcf --setFilteredGtToNocall -o FilteredStep2.vcf
java -Xmx8g -jar /opt/GenomeAnalysisTK/GenomeAnalysisTK.jar -T VariantsToTable -R Reference.fna -V FilteredStep2.vcf -F CHROM -F POS -GF GT -o FilteredFinalVCF.table

#Extract reference genotype.
cat -n pooledHC.vcf | awk '$2 == "#CHROM"' | awk '{print $1}' > StartNum.txt
for i in $(cat StartNum.txt); do tail -n +$((i)) pooledHC.vcf > CleanVCF.txt; done
awk '{print $4}' CleanVCF.txt > RefGT.txt
paste RefGT.txt FilteredFinalVCF.table > GtVcf.txt
rm StartNum.txt CleanVCF.txt RefGT.txt
#Extract data only for SNPs at the targeted 151st positions.
awk '$3 == "151"' GtVcf.txt > Pos151.txt
grep Target Pos151.txt > Pos151Target.txt
#Collect information for the SNPs.
sed -i 's/_Target//g' Pos151Target.txt
awk '{print $2}' Pos151Target.txt > List.txt
while read A; do grep -w "$A" 2770TargetInfo.txt >> ListInfo.txt; done < List.txt
cut -f 2,3 --complement Pos151Target.txt > Genotype.txt
wc -l ListInfo.txt Genotype.txt #Their numbers of lines should be the same.
paste ListInfo.txt Genotype.txt > a.txt
#Add headers and output results.
head -1 2770TargetInfo.txt > Header1.txt
head -1 GtVcf.txt | cut -f 2,3 --complement > Header2.txt
paste Header1.txt Header2.txt > Headers.txt
cat Headers.txt a.txt > Results.txt
rm Pos151.txt List.txt Pos151Target.txt ListInfo.txt Genotype.txt GtVcf.txt a.txt Header*

#The Results.txt should look like this:
#Target   Ref   Chr   Pos       REF   sample37.GT sample43.GT sample48.GT
#M_802    Tif1  A02   96307368   G       G/G         A/A         G/G
#M_2987   Tif1  A08   5056232    A       G/G         A/A          .
#M_7795   Tif1  B10   43563936   T       T/T         C/C         T/T
#RA_332   Tif2  A10   6375783    G       ./.         G/T         G/G
#RB_62    Tif2  B03   14730070   T       T/G         T/G         T/G
#Genotype "." means missing data, and "./." means low quality (filtered out).

#Download the Results.txt file to PC and open in excel.
#Replace the sample numbers with accession names for further analysis (for example: sample37=TxAG-6, sample43=43-09-03-02, and sample48=Tifrunner).