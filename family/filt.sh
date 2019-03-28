#!/usr/bin/env bash
set -ue
sample="$1"
ref=/nfs/brubeck.bx.psu.edu/scratch4/nick/refs/hg19/chr/chrM-rCRS.fa
~/code/indels/nvc-filter.py -r S -c 1000 -f 0.75 $sample/vars.vcf \
  > $sample/indels.vcf
~/code/indels/inspect-reads.py -tl -S $sample $sample/asm_filt.bam \
  -V $sample/indels.vcf -r $sample/asm.fa > $sample/indels.tsv
lastz $ref $sample/asm.fa > $sample/lift.lav
~/code/indels/quick-liftover.py $sample/lift.lav $sample/indels.tsv \
  > $sample/indels-lift.tsv
awk -F '\t' '$6 > 1000 && $8 > 1 && $20 <= 1 && $21 <= 1 {print}' \
  $sample/indels-lift.tsv > $sample/indels-lift-filt.tsv
