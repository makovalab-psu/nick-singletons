#!/usr/bin/env bash
set -ue
function fail {
  echo "$1" >&2
  exit 1
}
DunovoDir="$HOME/code/dunovo-clean"
Processes=30
SortMem=30G
TmpDir='../../tmp'
if [[ "$#" -lt 1 ]]; then
  fail "Usage: $(basename $0) sample_dir dist"
fi
sample_dir="$1"
dist="$2"
dir="$sample_dir/dunovo.corr/dist$dist"
if ! [[ -d "$dir" ]]; then
  mkdir -p "$dir"
fi
cd "$dir"

ref=$(ls ../../ref/*.fa | head -n 1)
if ! [[ -d "$TmpDir" ]]; then
  mkdir "$TmpDir"
fi

function paste_mates {
  fastq1="$1"
  fastq2="$2"
  if [[ $(basename "$fastq1" .gz) == $(basename "$fastq1") ]]; then
    paste "$fastq1" "$fastq2"
  else
    paste <(gunzip -c "$fastq1") <(gunzip -c "$fastq2")
  fi
}
# Find the input reads.
for name in renamed.full_1.fq renamed.full_1.fq.gz reads-renamed_1.fq.gz; do
  if [[ -s "../../$name" ]]; then
    fastq1="../../$name"
    fastq2=$(echo "$fastq1" | sed -E 's/_1\./_2./')
  fi
done
if ! [[ "$fastq1" ]]; then
  fail "Borked!"
fi
# Run Du Novo.
echo "Making families.."
paste_mates "$fastq1" "$fastq2" \
  | paste - - - - \
  | awk -f "$DunovoDir/make-barcodes.awk" \
  | sort -S "$SortMem" -T "$TmpDir" \
  > families.tsv
echo "Aligning barcodes.."
"$DunovoDir/baralign.sh" families.tsv refdir correct.sam
echo "Running corrected half of pipeline.."
"$DunovoDir/correct.py" -I --dist "$dist" --mapq 20 --pos 2 \
    families.tsv refdir/barcodes.fa correct.sam \
  | sort -S "$SortMem" -T "$TmpDir" \
  | tee -a families.corr.tsv \
  | "$DunovoDir/align-families.py" -I -p "$Processes" \
  | tee -a families.corr.msa.tsv \
  | "$DunovoDir/make-consensi.py" -r 3 -q 25 --cons-thres 0.7 --fastq-out 40 \
    --sscs1 sscs.corr.raw_1.fq --sscs2 sscs.corr.raw_2.fq \
    -1 duplex.corr.raw_1.fq -2 duplex.corr.raw_2.fq
echo "Trimming corrected output.."
~/code/makrutenko/trimmer.py --format fastq --filt-bases N --thres 0.3 --window 10 \
    --min-length 75 duplex.corr.raw_1.fq duplex.corr.raw_2.fq \
    duplex.corr.filt_1.fq duplex.corr.filt_2.fq
echo "Running uncorrected half of pipeline.."
"$DunovoDir/align-families.py" -I -p "$Processes" families.tsv \
  | tee -a families.uncorr.msa.tsv \
  | $DunovoDir/make-consensi.py -r 3 -q 25 --cons-thres 0.7 --fastq-out 40 \
    --sscs1 sscs.uncorr.raw_1.fq --sscs2 sscs.uncorr.raw_2.fq \
    -1 duplex.uncorr.raw_1.fq -2 duplex.uncorr.raw_2.fq
echo "Trimming uncorrected output.."
~/code/makrutenko/trimmer.py --format fastq --filt-bases N --thres 0.3 --window 10 \
    --min-length 75 duplex.uncorr.raw_1.fq duplex.uncorr.raw_2.fq \
    duplex.uncorr.filt_1.fq duplex.uncorr.filt_2.fq
# Align.
for corr in uncorr corr; do
  for filt in raw filt; do
    echo "Aligning $filt $corr.."
    bamtmp="duplex.$corr.$filt.tmp.bam"
    bwa mem -M -t 32 "$ref" "duplex.$corr.${filt}_1.fq" "duplex.$corr.${filt}_2.fq" \
      | samtools view -Sb - \
      > "$bamtmp"
    samtools sort -o "$bamtmp" dummy > "duplex.$corr.$filt.bam"
    rm "$bamtmp"
  done
done
