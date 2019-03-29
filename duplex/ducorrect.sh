cd $1
~/code/duplex/baralign.sh families.tsv barref barcodes.bam 24
samtools view barcodes.bam | ~/code/duplex/correct.py families.tsv > families.corrected.tsv
sort families.corrected.tsv > families.corrected.sorted.tsv
~/code/duplex/align_families.py -p 30 families.corrected.sorted.tsv \
  | tee -a families.msa.tsv \
  | ~/code/duplex/dunovo.py -r 3 -q 25 --sscs-file sscs.all.fa > duplex.fa
