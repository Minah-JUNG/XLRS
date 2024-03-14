#!/bin/bash

while read F; do

#1. Subset of VCF file

echo "Creating Subset VCF file"

bcftools view -S ./$F\.txt -v snps -O z -o XLRS_$F\.vcf.gz ./GenotypeGVCFs_XLRS.snp_filt_minDP10.annotated.vcf.gz

zcat XLRS_$F\.vcf.gz | sed 's/|//g' > input.vcf

echo "VCF stat calculation"

perl VCF_to_hapmap.pl input.vcf $F

#2. Filtering

echo "Filtering"

awk '($5==0 && $12 == 0 && $9 == 0 && $15 == "O"){print $0}' $F\_SNP_statistics.txt > $F\_snp_list.txt

done<./list.txt

