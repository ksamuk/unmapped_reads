#!/bin/bash

FQ=$1
BAM=$2

# print out mapped headers
echo "extracting unmapped headers from BAM file..."
samtools view -F 4 $BAM | awk '{ print ">"$1"" }' | sort -u > mapped.txt

# grap only the fq headers.
echo "building list of original FASTA headers..."
awk '(NR % 2 == 1)' $FQ | sed 's/[ ]\[.*\]//g' | sort -u  > reads.txt

# join and print only un-paired (-v)
echo "building list of unmapped FASTA headers..."

join -v 1 reads.txt mapped.txt > unmapped.txt
#grep -F -v -f mapped.txt reads.txt > unmapped.txt

# pull unmapped reads out of fasta
echo "writing unmapped reads to file..."
grep -F -f unmapped.txt --no-group-separator -A 1 batch_1.fa > unmapped_reads.fa

rm mapped.txt
rm reads.txt
rm unmapped.txt