#!/usr/bin/env bash
root="$1"
sample="$2"
fastq1=$root/data/filt/${sample}_1.fq
fastq2=$root/data/filt/${sample}_2.fq
echo $sample
set -x
mkdir -p $root/asm/$sample
mkdir -p $root/data/filt/$sample
java -jar ~/src/picard-tools-1.100/SamToFastq.jar INPUT=$root/data/map/$sample/filt.bam FASTQ=$fastq1 SECOND_END_FASTQ=$fastq2 VALIDATION_STRINGENCY=SILENT
spades.py -k 21,33,55,77,99,127 --careful -1 $fastq1 -2 $fastq2 -o $root/asm/$sample
