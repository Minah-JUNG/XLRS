#!/bin/bash
while read F; do

	bcftools filter -S . -i 'GT="1/1" && %QUAL>=30 && GQ>=30 && DP>=10' $F\_raw.vcf.gz | sed '/\.\/\./d' | gzip -c > $F\_bcf_filtered.vcf.gz &

done<./$1

	echo "running"  > message.txt
