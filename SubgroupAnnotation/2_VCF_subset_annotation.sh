#!/bin/bash

	 ls -rlht ./ | awk '{print $5 "\t" $9}' | grep -Ev "^0" | grep "txt" | cut -f 2 | sed 's/_/\t/g' | cut -f 1,2 > step2_work_list.txt

while IFS='	' read E F; do

	echo "$E $F is working"

	cut -f 1 ./$E\_$F\_snp_list.txt | sed 's/_/\t/g' > variant_list.txt

	bcftools view -T variant_list.txt -v snps -O z -o XLRS_$E\_$F\.vcf.gz ./XLRS_$E\.vcf.gz

	rm variant_list.txt

	zcat XLRS_$E\_$F\.vcf.gz | grep "#" > header.txt

	zcat XLRS_$E\_$F\.vcf.gz | grep -Ev "#" | sed 's/;ANN=/\t/g' | cut -f 1-8,10- > genotype.txt

	cat header.txt genotype.txt > input.vcf

	rm header.txt genotype.txt

	java -jar snpEff.jar HG38_Ensemble input.vcf | gzip -c > XLRS_$E\_$F\.vcf.annotated.gz

	mv snpEff_genes.txt $E\_$F\_snpEff_genes.txt

	perl Auto_filler_for_empty_position_in_lines.pl $E\_$F\_snpEff_genes.txt

	mv snpEff_summary.html $E\_$F\_snpEff_summary.html

	cut -f 1-8 Filled_$E\_$F\_snpEff_genes.txt | sed '1,2d' | sort -k2,2 -k 5,5nr -k 6,6nr -k 7,7nr -k 8,8nr -k3,3 > sorted_$E\_$F\_snpEff_genes.txt

	perl Get_top_line_for_sorted_column.pl sorted_$E\_$F\_snpEff_genes.txt

	sort -k 5,5nr -k6,6nr -k7,7nr -k 8,8nr Top_Selected_sorted_$E\_$F\_snpEff_genes.txt > Final_Top_Selected_sorted_$E\_$F\_snpEff_genes.txt

	rm XLRS_$E\_$F\.vcf.gz input.vcf

done<./step2_work_list.txt
