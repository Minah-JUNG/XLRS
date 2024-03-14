#!/bin/bash
while read F; do

	samtools sort -@ 10 -O bam -o sorted_$F\.bam ./$F\.bam

	java -jar -Xmx20g picard.jar MarkDuplicates INPUT=./sorted_$F\.bam OUTPUT=$F\_rmdup.bam REMOVE_DUPLICATES=true METRICS_FILE=$F\_metrics.txt AS=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

	java -jar -Xmx20g picard.jar FixMateInformation INPUT=./$F\_rmdup.bam OUTPUT=$F\_fxmt.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE

	java -jar -Xmx20g GenomeAnalysisTK.jar -T RealignerTargetCreator -R ./1_Reference/GCF_000001405.40_GRCh38.p14_genomic.fasta -I $F\_fxmt.bam -o $F\_realign.list

	java -jar -Xmx20g GenomeAnalysisTK.jar -T IndelRealigner -R ./1_Reference/GCF_000001405.40_GRCh38.p14_genomic.fasta -I ./$F\_fxmt.bam -targetIntervals $F\_realign.list -o $F\_realign.bam

	java -jar -Xmx20g picard.jar AddOrReplaceReadGroups INPUT=./$F\_realign.bam OUTPUT=./$F\_fixed.bam SORT_ORDER=coordinate RGID=$F RGLB=$F RGPL=illumina RGPU=$F RGSM=$F CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

	rm $F\_rmdup.bam $F\_metrics.txt $F\_fxmt.bam $F\_fxmt.bai $F\_realign.list $F\_realign.bam $F\_realign.bai 

done<./$1
