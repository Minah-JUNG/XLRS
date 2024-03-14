#!/bin/bash

while IFS='	' read F G; do

	echo "chr $F of $G is working"

	zcat ./GenotypeGVCFs_XLRS_chr$F\.$G\_filt.vcf.gz | sed 's/|/\//g' | sed 's/chr//g' > input.vcf

	bcftools view -m 2 -M 2 -O z -o Biallelic_input.vcf.gz ./input.vcf

	java -jar -Xmx10g snpEff.jar HG38_Ensemble Biallelic_input.vcf.gz | gzip -c > Biallelic_input.annotated.vcf.gz

	mv Biallelic_input.annotated.vcf.gz GenotypeGVCFs_XLRS_chr$F\.$G\_filt.annotated.vcf.gz

	mv snpEff_summary.html GenotypeGVCFs_XLRS_chr$F\.$G\_summary.html

	mv snpEff_genes.txt GenotypeGVCFs_XLRS_chr$F\.$G\_genes.txt

	mv GenotypeGVCFs_XLRS_chr$F\.$G\_summary.html GenotypeGVCFs_XLRS_chr$F\.$G\_genes.txt GenotypeGVCFs_XLRS_chr$F\.$G\_filt.annotated.vcf.gz Annotated_Variant_calling_result

done<./$1
