#!/bin/bash
# same as filter.sh, but this only filters for location.
# Note: requires ~/code/sam-uniq-bowtie2.awk
set -ue

OUTFILE=compiled-unfiltered.sam
FOLDERS='odt-total odt-brain tseq-total tseq-brain'

first=1
cd ~/scratch/maqc
for folder in $FOLDERS; do
  files=$(ls $folder/bam/*.bam)

  for file in $files; do
    if [ $first == 1 ]; then
      cmd="samtools view -H $file > $OUTFILE"
      echo '$ '$cmd; eval "$cmd"
      first=0
    fi
    cmd="samtools view -h $file chrM:2617-2617 >> $OUTFILE"
    echo '$ '$cmd; eval "$cmd"
  done
done
