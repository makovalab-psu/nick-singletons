set -ue
i=$1
function validate {
  dir=$1
  last_bar=$(cut -f 1 $dir/families.tsv | tail -n 1)
  if ! grep -q "$last_bar" $dir/refdir/barcodes.fa || \
     ! grep -q "$last_bar" $dir/refdir/barcodes-ref.fa; then
    echo "Error: Last barcode $last_bar missing from $dir/refdir." >&2
    return 1
  fi
  name=$(grep -B 1 $last_bar $dir/refdir/barcodes.fa | head -n 1 | tail -c +2)
  if ! [[ $(awkt '$1 == '$name $dir/correct.sam | wc -l) -gt 0 ]]; then
    echo "Last family $last_bar doesn't appear in the $dir/correct.sam alignment." >&2
    return 1
  fi
  return 0
}
if ! [[ -d "dunovo$i" ]]; then
  mkdir "dunovo$i"
fi
paste input$i/reads_1.fq input$i/reads_2.fq \
  | paste - - - - \
  | awk -f ~/code/dunovo-clean/make-barcodes.awk \
  | srun -C new -J sort1 -c 64 --mem 100G sort -S 80% --parallel 60 \
    -T /nfs/brubeck.bx.psu.edu/scratch5/nick/tmp \
  > dunovo$i/families.tsv
sleep 10
srun -C new --mem 100G ~/code/dunovo-clean/baralign.sh dunovo$i/families.tsv dunovo$i/refdir \
  dunovo$i/correct.sam
sleep 1m
validate dunovo$i
last_bar=$(cut -f 1 dunovo$i/families.tsv | tail -n 1)
srun -C new -J correct --mem 100G ~/code/dunovo-clean/correct.py --dist 3 --mapq 20 \
    --pos 2 dunovo$i/families.tsv dunovo$i/refdir/barcodes.fa dunovo$i/correct.sam \
  | srun -C new -J sort2 -c 64 --mem 100G sort -S 80% --parallel 60 \
    -T /nfs/brubeck.bx.psu.edu/scratch5/nick/tmp \
  | tee dunovo$i/families.corrected.tsv \
  | srun -C new -J align --mem 100G -c 64 ~/code/dunovo-clean/align-families.py \
    -p 63 -a kalign \
  | tee dunovo$i/families.msa.tsv \
  | srun -C new -J consensi --mem 100G -c 64 ~/code/dunovo-clean/make-consensi.py \
    -r 3 -q 0 --fastq-out 40 --sscs1 dunovo$i/sscs_1.fq --sscs2 dunovo$i/sscs_2.fq \
    -1 dunovo$i/duplex_1.fq -2 dunovo$i/duplex_2.fq
