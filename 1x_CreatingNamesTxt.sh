##!/bin/bash

##Author: Joy Sung

##If you don't have a file including the barcode information below, run commands in this session to create the required "Names.txt" file.
##Note, you'll need to modify your own information using the "vi" command.
##Before running these commands, create the "names.txt" file by following steps at the beginning of "1_Variant_Calling.sh" file.

##Extract RGPU (Read Group Platform Unit; eg. run barcodes) to make a new list; remember to modify this part based on your file names.
#while read A; do head -1 ../"$A"_L001_R1_001.fastq >> Heads.txt; done < names.txt
#awk '{print $2}' Heads.txt > Barcode.txt #Check whether modification is needed.
#sed -i 's/1:N:0://g' Barcode.txt #Check whether modification is needed.
#paste names.txt Barcode.txt > a.txt

##Use "vi" then type "i" to edit the file by adding the 2nd and the 3rd columns manually (see example below).
#vi a.txt
##You just need to fill in the 2nd and the 3rd columns, no need to sort them for now, we'll do it later.
##Eventually you'll want to have them sorted in the same order as the original sample numbers.
##Here's an example (tab delimited) with 4 columns: names, numbers, samples, RGPU (barcodes):
##LBK22-B01_S2     37     sample37     AATCCAGC
##LBK22-A01_S1     43     sample43     CGCTACAT
##LBK22-A06_S41    48     sample48     ACGCTTCT
##After editing the file, press "esc", then type ":wq!" to save (use ":q!" if you don't want to save it).

##Sort it by the sample numbers (the 2nd column).
#sort -k2n a.txt > Names.txt
#wc Names.txt #This should print the number of samples.
#rm Heads.txt Barcode.txt a.txt

##Make sure all the information and format are correct in the "Names.txt" before proceeding to the next step.