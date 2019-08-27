#!/bin/bash
#
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20000

#############
## Read Me ##
#############
# It takes as input a bam file and the name of the dir where output files will be created.
# It will convert a bam file to fastq file, will discard reads shorter than `MinLen`.
#
# WARNING: Pipeline assumes SOLEXA QUALS!!! (if phred+64 quals remove -Q33 parameter got fastx_clipper)
#
# Last reviewed: 03/jul/2019	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      Bam2Fq -i BAM_FILE -o OUT_DIR [-l MIN_LENGTH]"
    exit 1
}
# Defaults
MinLen=20
# Parse options.
while getopts 'i:o:l:' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;				# BAM file
    o) outDir="${OPTARG}" ;;		# output directory.
    l) MinLen="${OPTARG}" ;;		# minimum read length (20)
    *) usage ;;
  esac
done

# Clip path and suffixes from filename.
prefix=$(basename $bam)
prefix=${prefix/.bam/''}
prefix=${prefix/.sorted/''}
prefix=$(perl -e 'if($ARGV[0]=~/(.+)_\d{8}\w?_\d{8}$/){print $1}else{print $ARGV[0]}' $prefix)   # Get rid of the date, to make names shorter and keep them manageable when more things are appended to them.
####################


## Module Loading ##
module load fastx-toolkit/0.0.14-foss-2017a
module load samtools/1.9-foss-2017a
#####################

# Create destination.
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
  mkdir -p ${outDir}
fi

cleanup=0
# If BAM not sorted, sort it.
samtools view -H $bam | head -n1 > ${outDir}/${prefix}.temp_header.txt
if grep -E -q "SO:coordinate" ${outDir}/${prefix}.temp_header.txt; then
	echo "${prefix}: sorted bam as input"
else
	echo "${prefix}: sorting input bam"
	samtools sort -@ 3 -o ${outDir}/${prefix}.sorted.bam $bam
	bam="${outDir}/${prefix}.sorted.bam"
  cleanup=1
fi
rm ${outDir}/${prefix}.temp_header.txt

# BAM to FASTQ.
echo "${prefix}: Extracting reads"
samtools bam2fq -@ 3 -0 ${outDir}/${prefix}.MyFq.fq ${bam}
# Remove reads that are too short.
echo "${prefix}: Filtering length"
fastx_clipper -d200 -Q33 -l${MinLen} -i ${outDir}/${prefix}.MyFq.fq -o ${outDir}/${prefix}.fq
rm ${outDir}/${prefix}.MyFq.fq
echo "${prefix}: Compressing file"
gzip -f ${outDir}/${prefix}.fq

# Cleanup
if [ $cleanup -eq 1 ]; then
  rm $bam
fi

echo "${prefix}: Finished read extraction and length filtering"
exit $?
####################
