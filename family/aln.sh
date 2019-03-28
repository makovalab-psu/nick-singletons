#!/usr/bin/env bash
set -ue
REF=/nfs/brubeck.bx.psu.edu/scratch2/nick/refs/spike/bwa/hg19-rCRS-phiX-pUC18.fa
cd /nfs/brubeck.bx.psu.edu/scratch2/nick/family2/20141001-TR3
fastq1=$1
if ! echo $fastq1 | grep -q _R1_ ; then
  exit
fi
fastq2=$(echo $fastq1 | sed 's/_R1_/_R2_/')
base=$(basename $fastq1 .fastq.gz | sed 's/_[12]//' | sed -E 's/_S[0-9]+_L001_R[12]_001//')
echo $base
bwa mem -t 16 -M $REF data/$fastq1 data/$fastq2 > aln/$base.sam
samtools view -Sb aln/$base.sam > aln/$base.tmp.bam
samtools sort aln/$base.tmp.bam aln/$base
samtools index aln/$base.bam
rm aln/$base.sam aln/$base.tmp.bam
