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
#SBATCH --qos=medium
#SBATCH --time=1-0:00:00       # one day

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i BAM -r REF_FASTA -b UTR_BED -a GTF -o OUT_DIR -q MINBQ -H HARD_CLIP"
    exit 1
}
# Parse options.
minq=10
adaptor='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'  # read 1, single end 5'adaptor
while getopts 'i:r:b:a:o:q:H:' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;			# unaligned bam (one file)
    r) ref="${OPTARG}" ;;			# whole reference fasta for alignment
    b) bed="${OPTARG}" ;;			# bed of 3'UTRs (for map filtering)
    a) annot="${OPTARG}" ;;			# GTF annotation (for counting)
    o) outdir="${OPTARG}" ;;		# output directory
    q) minq="${OPTARG}" ;;			# minimum alignement quality (10)
    H) trim3="${OPTARG}" ;;         # trim this many bases from end of reads
    *) usage ;;
  esac
done

if [ ! -z "$trim3" ]; then
    trim3="-u -${trim3}"
fi

prefix=$(basename $bam)
prefix=${prefix/.bam/}

mkdir -p ${outdir}/fastqc_pre ${outdir}/fastqc_post ${outdir}/fastq ${outdir}/fastq_trimmed ${outdir}/mapped ${outdir}/mapped_filtered ${outdir}/counts ${outdir}/multiqc


source ~/miniconda3/bin/deactivate
source ~/miniconda3/bin/activate bioinfo # For me it's always loaded


echo "$bam - BAM to FASTQ"
samtools bam2fq -0 ${outdir}/fastq/${prefix}.fastq $bam
gzip -f ${outdir}/fastq/${prefix}.fastq

echo "$bam - QC"
fastqc -o ${outdir}/fastqc_pre ${outdir}/fastq/${prefix}.fastq.gz

echo "$bam - trim reads and filter remaining lengths"
cutadapt -j 4 -a A{18} -a T{18} -g $adaptor -u +12 ${trim3} -q 20 -m 20 -o ${outdir}/fastq_trimmed/${prefix}.fastq.gz ${outdir}/fastq/${prefix}.fastq.gz > ${outdir}/fastq_trimmed/${prefix}.cutadapt.log

echo "$bam - QC"
fastqc -o ${outdir}/fastqc_post ${outdir}/fastq_trimmed/${prefix}.fastq.gz


source ~/miniconda3/bin/deactivate
source ~/miniconda3/bin/activate slamdunk


echo "$bam - map"
slamdunk map -t 4 -r $ref -o ${outdir}/mapped -5 0 -n 100 --quantseq -ss ${outdir}/fastq_trimmed/${prefix}.fastq.gz
echo "$bam - filter locations"
slamdunk filter -t 4 -o ${outdir}/mapped_filtered -b $bed ${outdir}/mapped/${prefix}.fastq_slamdunk_mapped.bam


source ~/miniconda3/bin/deactivate
source ~/miniconda3/bin/activate bioinfo


echo "$bam - quantify alignments with qual >$minq"
mkdir -p ${outdir}/counts
featureCounts -T 4 -M --primary -Q $minq -a $annot -o ${outdir}/counts/${prefix}.txt ${outdir}/mapped_filtered/${prefix}.fastq_slamdunk_mapped_filtered.bam


#
# The steps below don't make sense for individual files. They work best after the previous steps are all finished for all files.
#

# echo "$bam - MultiQC"
# srun multiqc -f -o ${outdir}/multiqc ${outdir}/fastq ${outdir}/fastq_trimmed ${outdir}/fastqc_pre ${outdir}/fastqc_post ${outdir}/mapped ${outdir}/mapped_filtered ${outdir}/counts

# Extract and collate the counts, lose the commented command.
# fileutilities.py T ${outdir}/counts/*txt -i -l --cols 6 | perl -e 'while(<>){~s/_\|6//g;print}' > ${outdir}/counts/all_counts.tsv
