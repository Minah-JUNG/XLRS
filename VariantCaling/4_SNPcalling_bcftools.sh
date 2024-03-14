#!/bin/bash
while read F; do

	bcftools mpileup --threads 40 -O u -f ./1_Reference/GCF_000001405.40_GRCh38.p14_genomic.fasta ./3_Fixed_bam/$F\_fixed.bam | bcftools call --threads 40 -f GQ -V indels -vmO z > $F\_raw.vcf.gz &

done<./$1

	echo "work"
