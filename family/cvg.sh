#!/usr/bin/env bash
bamfile=$1
outfile=$2
base=$(basename $bamfile .bam | cut -f 1 -d _ | cut -f 1 -d -)
echo $base
samtools view -b $bamfile chrM | bamtools coverage | awk -v OFS='\t' '{print "'$base'", $2, $3}' >> $outfile

