#!/bin/bash
while read F; do

	mkdir $F

        java -jar -Xmx100g trimmomatic-0.39.jar PE -phred33 -threads 40 ./Input_symbolic/$F\_1.fq.gz ./Input_symbolic/$F\_2.fq.gz ./$F/$F\_1.fq.gz ./$F/Unpaired_$F\_1.fq.gz ./$F/$F\_2.fq.gz ./$F/Unpaired_$F\_2.fq.gz  ILLUMINACLIP:/home/sungmin716/Programs/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:75

done<./$1
