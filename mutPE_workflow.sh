#!/usr/bin/env sh

# Process the samples of one Run, from .sra to _barcode_report.pdf

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -d BASEDIR -b BATCHDIR -i BOWTIE2_IDX [-D DATADIR] [-p PROCESSING_DIR] [-r RESULTS_DIR] [-a AUX_DIR] [[-s] | [-f]]"
		exit 1
}
# Defaults
data='data'
process='process'
results='results'
aux='aux'
# Parse options.
while getopts 'd:D:b:p:r:a:i:sf' flag; do
  case "${flag}" in
    d) base="${OPTARG}" ;;        # Base dir
		D) data="${OPTARG}" ;;     		# Data dir in which the batch is located, relative to base
		b) run="${OPTARG}" ;;         # Batch dir, relative to data dir
		p) process="${OPTARG}" ;;     # Dir in which intermediate steps output files are created, relative to base
		r) results="${OPTARG}" ;;     # Dir for the final files, relative to base
		a) aux="${OPTARG}" ;;     		# Dir where the bowtie index is located, relative to base
    i) bowtie2idx="${OPTARG}" ;;  # Bowtie2 index prefix
		s) issra="${OPTARG}" ;;       # Input is .sra format
		f) renamefq="${OPTARG}" ;;    # Files are uncompressed .fq instead of compressed .fastq.gz
		*) usage ;;
  esac
done


echo "Creating destinations."
if [[ ! -d "${base}/${process}/${run}" ]]; then
	mkdir -p "${base}/${process}/${run}"
fi
if [[ ! -d "${base}/${results}/${run}" ]]; then
	mkdir -p "${base}/${results}/${run}"
fi


if [[ -z $renamefq ]]; then
	echo "Renaming & compressing .fq to .fastq.gz."
	fileutilities.py T ${base}/${data}/${run} --dir fq | fileutilities.py P --loop S mv {abs} {dir}/{bas}.fastq \; srun gzip {dir}/{bas}.fastq
fi

if [[ -z $issra ]]; then
	echo "Converting .sra to fastq.gz."
	fileutilities.py T ${base}/${data}/${run} --dir .sra | fileutilities.py P --loop S srun ~/utility_scripts/sra2fastq.sh ${base}/${data}/${run} ${base}/${data}/${run}/{val}
fi

echo "Merging pairs."
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop S srun flash ,-t 1 ,-M 150 ,-x 0.25 ,-d ${base}/${process}/${run} ,-o {ali}_x25 ,-z ${base}/${data}/${run}/{ali}_1.fastq.gz ${base}/${data}/${run}/{ali}_2.fastq.gz

echo "Aligning."
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop S srun ,-c 4 ~/utility_scripts/bowtie2aln.sh ${base}/${aux}/${bowtie2idx} ${base}/${process}/${run} ${base}/${process}/${run} {ali}_x25.extendedFrags false 4 ,--very-sensitive-local


echo "Stratifying BAMs by number of mismatches."
fileutilities.py T ${base}/${process}/${run} --dir 'sorted.bam$' | fileutilities.py P --loop S srun ~/utility_scripts/split_BAM_by_numMismatches.sh {abs} 1 2
fileutilities.py T ${base}/${process}/${run} --dir 'sorted.bam$' | fileutilities.py P --loop S srun ~/utility_scripts/split_BAM_by_numMismatches.sh {abs} 3 10
fileutilities.py T ${base}/${process}/${run} --dir 'sorted.bam$' | fileutilities.py P --loop S srun ~/utility_scripts/split_BAM_by_numMismatches.sh {abs} 11 30


echo "Piling up."
fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.bam$' | fileutilities.py P --loop S srun samtools mpileup ,-A ,-B ,-f ${base}/${aux}/${bowtie2idx}.fa ,-d 1000000000 {abs} \> {dir}/{bas}.pileup

echo "Counting."
fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.pileup$' | fileutilities.py P --loop S python3 ./analysis/MutPE_quantification.py ,-p {abs} \> {dir}/{bas}.stats

echo "Visualising."
mutPE_mutation-stats_viz.R ${base}/${results}/${run} ${run} NULL yes ${base}/${process}/${run}/*.stats

echo "Finished."
