#!/bin/bash

# Define the base directories and references for easier access and modification
DirInput="./Raw"
DirRef="./Reference_Ensembl/GRCh38"
Ref="$DirRef/Homo_sapiens.GRCh38.dna.chromosome.fa"
DirFastq="./1_Trimmed_fastq"
DirBam="./2_BamSam"
DirFixedBam="./3_Fixed_bam"
DirRawVcf="./4_Raw_VCF"
DirFilteredVcf="./5_Filtered_VCF"
DirAnnotatedVcf="./6_Annotated_VCF" 

## making directories
mkdir -p "$DirFastq" "$DirBam" "$DirFixedBam" "$DirRawVcf" "$DirFilteredVcf" "$DirAnnotatedVcf"

## Variant Calling
while read -r F; do
    ## Step 1: Read trimming with Trimmomatic
    java -jar -Xmx100g trimmomatic-0.39.jar PE -phred33 -threads 40 "$DirInput/$F_1.fq.gz" "$DirInput/$F_2.fq.gz" "$DirFastq/$F_1.trimmed.fq.gz" "$DirFastq/Unpaired_$F_1.fq.gz" "$DirFastq/$F_2.trimmed.fq.gz" "$DirFastq/Unpaired_$F_2.fq.gz" ILLUMINACLIP:/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:75

    ## Step 2: Align reads and convert SAM to BAM
    bwa mem -M -R "@RG\tID:$F\tLB:$F\tPL:ILLUMINA\tSM:$F" -t 20 "$Ref" "$DirFastq/$F_1.trimmed.fq.gz" "$DirFastq/$F_2.trimmed.fq.gz" > "$DirBam/$F.sam"
    samtools view -@ 40 -b -q 10 -F 4 "$DirBam/$F.sam" > "$DirBam/$F.bam"

    ## Step 3: Sort BAM and mark duplicates
    samtools sort -@ 10 -O bam -o "$DirBam/sorted_$F.bam" "$DirBam/$F.bam"
    java -jar -Xmx20g picard.jar MarkDuplicates INPUT="$DirBam/sorted_$F.bam" OUTPUT="$DirFixedBam/$F_rmdup.bam" METRICS_FILE="$DirFixedBam/$F_metrics.txt" REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
    java -jar -Xmx20g picard.jar FixMateInformation INPUT="$DirFixedBam/$F_rmdup.bam" OUTPUT="$DirFixedBam/$F_fxmt.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE

    ## Step 4: Realigner and read group adjustment
    java -jar -Xmx20g GenomeAnalysisTK.jar -T RealignerTargetCreator -R "$Ref" -I "$DirFixedBam/$F_fxmt.bam" -o "$DirFixedBam/$F_realign.list"
    java -jar -Xmx20g GenomeAnalysisTK.jar -T IndelRealigner -R "$Ref" -I "$DirFixedBam/$F_fxmt.bam" -targetIntervals "$DirFixedBam/$F_realign.list" -o "$DirFixedBam/$F_realign.bam"
	java -jar -Xmx20g picard.jar AddOrReplaceReadGroups INPUT="$DirFixedBam/$F_realign.bam" OUTPUT="$DirFixedBam/$F_fixed.bam" SORT_ORDER=coordinate RGID="$F" RGLB="$F" RGPL=illumina RGPU="$F" RGSM="$F" CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

	## Step 5: Variant calling and filtering
	bcftools mpileup --threads 40 -Ou -f "$Ref" "$DirFixedBam/$F_fixed.bam" | bcftools call --threads 40 -mO z -o "$DirRawVcf/$F_raw.vcf.gz"
	bcftools filter -S . -i 'GT="1/1" && %QUAL>=30 && GQ>=30 && DP>=10' "$DirRawVcf/$F_raw.vcf.gz" | sed '/\.\/\./d' | gzip -c > "$DirFilteredVcf/$F_bcf_filtered.vcf.gz"

	## Step 6: Annotate variants
	java -jar -Xmx10g snpEff.jar HG38_Ensemble "$DirFilteredVcf/$F_bcf_filtered.vcf.gz" | gzip -c > "$DirAnnotatedVcf/$F_bcf_filtered.annotated.vcf.gz"

	#rm -f "$DirBam/$F.sam" "$DirBam/$F.bam" "$DirBam/sorted_$F.bam" "$DirFixedBam/$F_rmdup.bam" "$DirFixedBam/$F_fxmt.bam" "$DirFixedBam/$F_realign.bam"

	echo "$F processing complete"
	
done < "./$1"

##
echo "All processing finished."