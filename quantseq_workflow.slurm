#!/bin/bash
#
#SBATCH --get-user-env
#SBATCH -J qntSeq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20000
#SBATCH --output=qntSeq.out
#SBATCH --error=qntSeq.err

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i BAM -r REF_FASTA -b UTR_BED -a GTF -o OUT_DIR [-f FEATURE]"
    exit 1
}
# Parse options.
feature='exon'
minq=10
while getopts 'i:r:b:a:o:f:q:' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;			# unaligned bam (one file)
    r) ref="${OPTARG}" ;;			# whole reference fasta for alignment
    b) bed="${OPTARG}" ;;			# bed of 3'UTRs (for map filtering)
    a) annot="${OPTARG}" ;;			# GTF annotation (for counting)
		o) outdir="${OPTARG}" ;;			# output directory
    f) feature="${OPTARG}" ;;			# feature to quantify (gene)
    q) minq="${OPTARG}" ;;			# minimum alignement quality (0)
    *) usage ;;
  esac
done

prefix=$(basename $bam)
prefix=${prefix/.bam/}

source ~/miniconda3/bin/activate py27
module load fastqc/0.11.5-java-1.8.0_121 

echo $bam
echo "$bam - BAM to FASTQ"
samtools bam2fq -0 ${outdir}/fastq/${prefix}.fastq
gzip ${outdir}/fastq/${prefix}.fastq

echo "$bam - trim and FastQC"
trim_galore --fastqc --gzip --stringency 3 -o ${outdir}/fastq_trimmed ${outdir}/fastq/${prefix}.fastq.gz

echo "$bam - align"
slamdunk map -r $ref -o ${outdir}/mapped -5 12 -n 100 -t 4 --quantseq -ss ${outdir}/fastq_trimmed/${prefix}_trimmed.fq.gz
echo "$bam - filter"
slamdunk filter -o ${outdir}/map_filtered -b $bed -t 4 ${outdir}/mapped/${prefix}_trimmed.fq_slamdunk_mapped.bam

conda deactivate
source ~/miniconda3/bin/activate mybasics

echo "$bam - quantify alignments with qual >$minq"
mkdir -p ${outdir}/counts
htseq-count -f bam -s yes --nonunique all -t $feature -a $minq ${outdir}/map_filtered/${prefix}_trimmed.fq_slamdunk_mapped_filtered.bam $annot > ${outdir}/counts/${prefix}.txt

conda deactivate
module unload fastqc
