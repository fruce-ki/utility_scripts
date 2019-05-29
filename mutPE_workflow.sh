#!/usr/bin/env sh

# Process the samples of one Run, from .sra to _barcode_report.pdf

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -d BASEDIR -b BATCHDIR -i BOWTIE2_IDX [-D DATADIR] [-p PROCESSING_DIR] [-r RESULTS_DIR] [-a AUX_DIR] [-o OFFSET] [[-s] | [-f]]"
		exit 1
}
# Defaults
data='data'
process='process'
results='results'
aux='aux'
offsets='-:0:0'
# Parse options.
while getopts 'd:D:b:p:r:a:i:o:l:sf' flag; do
  case "${flag}" in
    d) base="${OPTARG}" ;;        # Base dir
		D) data="${OPTARG}" ;;     		# Data dir in which the batch is located, relative to base
		b) run="${OPTARG}" ;;         # Batch dir, relative to data dir
		p) process="${OPTARG}" ;;     # Dir in which intermediate steps output files are created, relative to base
		r) results="${OPTARG}" ;;     # Dir for the final files, relative to base
		a) aux="${OPTARG}" ;;     		# Dir where the bowtie index is located, relative to base
    i) bowtie2idx="${OPTARG}" ;;  # Bowtie2 index prefix
    o) offsets="${OPTARG}" ;;     # Reference deletions: "REF:START:LENGTH" ie. "HDR2:280:6,HDR2:280:6"
    s) issra="${OPTARG}" ;;       # Input is .sra format
		f) renamefq="${OPTARG}" ;;    # Files are uncompressed .fq instead of compressed .fastq.gz
		*) usage ;;
  esac
done


module load pandoc/2.0.5
module load flash/1.2.11-foss-2017a
module load samtools/1.9-foss-2017a


if [[ ! -d "${base}/${process}/${run}" ]]; then
	mkdir -p "${base}/${process}/${run}"
  echo "Created destination: ${base}/${process}/${run}"
fi
if [[ ! -d "${base}/${results}/${run}" ]]; then
	mkdir -p "${base}/${results}/${run}"
  echo "Created destination: ${base}/${results}/${run}"
fi


if [[ $renamefq ]]; then
	echo "Renaming & compressing .fq to .fastq.gz"
	fileutilities.py T ${base}/${data}/${run} --dir fq | fileutilities.py P --loop mv {abs} {dir}/{bas}.fastq \; srun gzip {dir}/{bas}.fastq
fi

if [[ $issra ]]; then
	echo "Converting ${base}/${data}/${run}/*.sra to ${base}/${data}/${run}/*.fastq.gz"
	fileutilities.py T ${base}/${data}/${run} --dir .sra | fileutilities.py P --loop srun ~/utility_scripts/sra2fastq.sh ${base}/${data}/${run} ${base}/${data}/${run}/{val}
fi

echo "Extracting ${base}/${data}/${run}/*.bam to ${base}/${data}/${run}/*_1/2.fastq.gz"
fileutilities.py T ${base}/${data}/${run} --dir .bam | fileutilities.py P --loop srun samtools fastq ,-c 1 ,-1 {dir}/{ali}_1.fastq.gz ,-2 {dir}/{ali}_2.fastq.gz {abs}

echo "Merging pairs in ${base}/${data}/${run} to ${base}/${process}/${run}/*_x25..fastq.gz"
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop srun flash ,-t 1 ,-M 150 ,-x 0.25 ,-d ${base}/${process}/${run} ,-o {ali}_x25 ,-z ${base}/${data}/${run}/{ali}_1.fastq.gz ${base}/${data}/${run}/{ali}_2.fastq.gz

echo "Aligning ${base}/${data}/${run}/*.fastq.gz to ${base}/${aux}/${bowtie2idx} and outputting to ${base}/${process}/${run}/*_x25.extendedFrags.bam"
fileutilities.py T ${base}/${data}/${run} --dir _1.fastq.gz | perl -e 'while(<>){~s/_1\.fastq$//;print}' | fileutilities.py P --loop srun ,-c 4 ~/utility_scripts/bowtie2aln.sh ${base}/${aux}/${bowtie2idx} ${base}/${process}/${run} ${base}/${process}/${run} {ali}_x25.extendedFrags false 4 ,--very-sensitive-local


echo "Stratifying ${base}/${process}/${run}/*.sorted.bam by number of mismatches"
fileutilities.py T ${base}/${process}/${run} --dir 'sorted.bam$' | fileutilities.py P --loop srun ~/utility_scripts/split_BAM_by_numMismatches.sh {abs} 1 2
fileutilities.py T ${base}/${process}/${run} --dir 'sorted.bam$' | fileutilities.py P --loop srun ~/utility_scripts/split_BAM_by_numMismatches.sh {abs} 3 10
fileutilities.py T ${base}/${process}/${run} --dir 'sorted.bam$' | fileutilities.py P --loop srun ~/utility_scripts/split_BAM_by_numMismatches.sh {abs} 11 30


echo "Piling up ${base}/${process}/${run}/\d+-\d+.bam onto ${base}/${aux}/${bowtie2idx}.fa"
fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.bam$' | fileutilities.py P --loop srun samtools mpileup ,-A ,-B ,-f ${base}/${aux}/${bowtie2idx}.fa ,-d 1000000000 {abs} \> {dir}/{bas}.pileup

echo "Counting."
fileutilities.py T ${base}/${process}/${run} --dir '\d+-\d+.pileup$' | fileutilities.py P --loop python3 ./analysis/MutPE_quantification.py ,-p {abs} \> {dir}/{bas}.stats

echo "Visualising ${base}/${process}/${run}/*.stats into ${base}/${results}/${run}/${run//\//_}*.html/pdf"
mutPE_mutation-stats_viz.R ${base}/${results}/${run} ${run/\//_} NULL no $offsets ${base}/${process}/${run}/*.stats

echo "Finished."
