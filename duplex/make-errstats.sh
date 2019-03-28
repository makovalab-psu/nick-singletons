set -ue
families="$1"
sscs_base="$2"
ref="$3"
outdir="$4"
input_dir=$(dirname "$sscs_base")
# Align.
sscs1="${sscs_base}_1.fq"
sscs2="${sscs_base}_2.fq"
if ! [[ -s "$sscs1" ]] || ! [[ -f "$sscs1" ]]; then
  echo "Error: sscs fastq file(s) missing: ${sscs_base}_[12].fq" >&2
  exit 1
fi
bwa mem -M -t 16 "$ref" "$sscs1" "$sscs2" > "$input_dir/sscs.sam"
samtools view -Sb "$input_dir/sscs.sam" > "$input_dir/sscs.tmp.bam"
samtools sort "$input_dir/sscs.tmp.bam" "$input_dir/sscs"
rm "$input_dir/sscs.tmp.bam" "$input_dir/sscs.sam"
# Run errstats.py.
mkdir -p "$outdir"
PYTHONPATH=~/code/pyBamParser/lib ~/code/dunovo-clean/utils/errstats.py -v --out-format errors2 \
    --no-indels --dedup --qual-errors -q 25 --min-reads 5 \
    --overlap-stats "$outdir/overlaps.snvs.tsv" --dedup-log "$outdir/dedup.snvs.log" \
    --bam "$input_dir/sscs.bam" "$families" \
  | gzip -c - \
  > "$outdir/errstats.snvs.tsv.gz"
# Collate errstats.
gunzip -c "$outdir/errstats.snvs.tsv.gz" \
  | gawk -F '\t' -v OFS='\t' '{for (i=8; i<=NF; i++) {tot[$4][$i]++}}
          END {for (f in tot) {for (e in tot[f]) {print f, e, tot[f][e]}}}' \
  > "$outdir/errstats.snvs.supercollated.tsv"
