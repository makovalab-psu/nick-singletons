#!/bin/bash
# Note: requires ~/code/sam-uniq-bowtie2.awk
set -ue

folders='odt-total odt-brain tseq-total tseq-brain'

for folder in $folders; do
  cd ~/scratch/maqc/$folder
  files=$(ls bam/*.bam)
  if [ ! -d filt ]; then
    mkdir filt
  fi

  for file in $files; do
    base=$(echo $file | sed -r 's/^.*\/([^/]*)\.bam$/\1/')
    cmd="samtools view -h -F 4 -q 1 $file chrM:2617-2617 | awk -f ~/code/sam-uniq-bowtie2.awk > filt/$base.sam"
    echo '$ '$cmd; eval "$cmd"
  done
done
