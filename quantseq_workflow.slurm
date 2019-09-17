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
    echo "      $0 -i BAM -r REF_FASTA -b UTR_BED -a GTF -o OUT_DIR"
    exit 1
}
# Parse options.
minq=10
while getopts 'i:r:b:a:o:q:' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;			# unaligned bam (one file)
    r) ref="${OPTARG}" ;;			# whole reference fasta for alignment
    b) bed="${OPTARG}" ;;			# bed of 3'UTRs (for map filtering)
    a) annot="${OPTARG}" ;;			# GTF annotation (for counting)
    o) outdir="${OPTARG}" ;;			# output directory
    q) minq="${OPTARG}" ;;			# minimum alignement quality (10)
    *) usage ;;
  esac
done

prefix=$(basename $bam)
prefix=${prefix/.bam/}

source ~/miniconda3/bin/activate py27
module load fastqc/0.11.5-java-1.8.0_121
module load subread/1.6.4-foss-2017a 

mkdir -p ${outdir}/fastq ${outdir}/mapped ${outdir}/map_filtered ${outdir}/counts ${outdir}/multiqc

echo "$bam - BAM to FASTQ"
samtools bam2fq -0 ${outdir}/fastq/${prefix}.fastq $bam
gzip ${outdir}/fastq/${prefix}.fastq

echo "$bam - trim and FastQC"
trim_galore --fastqc --gzip --stringency 3 -o ${outdir}/fastq_trimmed ${outdir}/fastq/${prefix}.fastq.gz

echo "$bam - align"
slamdunk map -t 4 -r $ref -o ${outdir}/mapped -5 12 -n 100 --quantseq -ss ${outdir}/fastq_trimmed/${prefix}_trimmed.fq.gz
echo "$bam - filter"
slamdunk filter -t 4 -o ${outdir}/map_filtered -b $bed ${outdir}/mapped/${prefix}_trimmed.fq_slamdunk_mapped.bam

echo "$bam - quantify alignments with qual >$minq"
mkdir -p ${outdir}/counts
featureCounts -T 4 -M --primary -Q $minq -a $annot -o ${outdir}/counts/${prefix}.txt ${outdir}/map_filtered/${prefix}_trimmed.fq_slamdunk_mapped_filtered.bam

conda deactivate
module unload fastqc
module unload subread

echo "$bam - MultiQC"
module load multiqc/1.3-foss-2017a-python-2.7.13
srun multiqc -f -o ${outdir}/multiqc ${outdir}/fastq ${outdir}/fastq_trimmed ${outdir}/mapped ${outdir}/map_filtered ${outdir}/counts
module unload multiqc

source ~/miniconda3/bin/activate mybasics
# Extract and collate the counts, lose the commented command.
#fileutilities.py T ${outdir}/counts/*txt --loop cut ,-f 1,7 {abs} \> {dir}/{bas}.tsv
fileutilities.py T ${outdir}/counts/*txt --loop cut ,-f 1,7 {abs} \| tail ,-n+2 \> {dir}/{bas}.tsv
fileutilities.py T ${outdir}/counts/*tsv -i -l --appnd outer > ${outdir}/counts/all_counts.txt
head -n2 ./process/counts/all_counts.txt | perl -e 'while(<>){~s/CD4LUANXX_4_//g;~s/_\w{6}_\|1//g;print}' > tmp ; tail -n +2 ./process/counts/all_counts.tsv >> tmp ; mv tmp ./process/counts/all_counts.tsv
