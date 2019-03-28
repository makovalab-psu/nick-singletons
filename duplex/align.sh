set -ue
dir=$1
ref=$2
bwa mem -M -t 16 $ref $dir/sscs_1.fq $dir/sscs_2.fq > $dir/sscs.sam
samtools view -Sb $dir/sscs.sam > $dir/sscs.tmp.bam
samtools sort $dir/sscs.tmp.bam $dir/sscs
rm $dir/sscs.tmp.bam $dir/sscs.sam
