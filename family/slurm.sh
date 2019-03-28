base=$1
cd /nfs/brubeck.bx.psu.edu/scratch2/nick
export LC_ALL=C
mkdir -p asm/spades/$base/merge2
mkdir -p asm/spades/$base/fastq
fqpair2fasta.pl family/all1/data/fastq/${base}_1.fq family/all1/data/fastq/${base}_2.fq asm/spades/$base/fastq/$base
contigMerger.pl asm/spades/$base/orient/contigs_orientedContigs hg19/chr/chrM-rCRS.fa asm/spades/$base/merge2/contig -mincontlen 200 -readfa asm/spades/$base/fastq/$base.fa -readq asm/spades/$base/fastq/$base.qual
