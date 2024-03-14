#!/bin/bash
while read F; do

	bwa mem -M -R "@RG\tID:ID\tLB:LB\tPL:ILLUMINA\tSM:SM" -t 20 ./1_Reference_Ensembl/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.fa ./2_Fastq/$F/$F\_1.fq.gz ./2_Fastq/$F/$F\_2.fq.gz > $F\.sam

	samtools view -@ 40 -b -q 10 -F 4 $F\.sam > $F\.bam

done<./$1
