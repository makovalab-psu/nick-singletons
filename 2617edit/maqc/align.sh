#!/bin/bash
set -ue

cd ~/scratch/maqc
folders='odt-total odt-brain tseq-total tseq-brain'

for folder in $folders; do
  files=$(ls $folder/fastq)

  for file in $files; do
    file=$(echo $file | sed -r 's/\.fastq$//')
    cmd="bowtie2 -p 3 --trim5 4 -x ~/scratch/hg19/whole/hg19 -U $folder/fastq/$file.fastq > $folder/sam/$file.sam"
    echo '$ '$cmd; eval "$cmd"
  done
done
