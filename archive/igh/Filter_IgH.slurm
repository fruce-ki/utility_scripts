#!/bin/bash
#
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40000

#############
## Read Me ##
#############
# It takes as input a (deduplicated) aligned BAM and the name of the dir where output files will be created.
#
# Last reviewed: 05/aug/2019	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -b BAM_FILE -o OUTDIR -c CELLLINE [-x SCRIPTDIR]"
    exit 1
}
# Defaults.
SCRIPT_PATH='/groups/pavri/Kimon/ighrealignment'
# Parse options.
while getopts 'b:o:c:l:x:' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;           # BAM file.
    o) outDir="${OPTARG}" ;;        # Output directory.
    c) CELLLINE="${OPTARG}" ;;      # Full reference name of locus chromosome.
    x) SCRIPT_PATH="${OPTARG}" ;;   # Directory with the python scripts (ie. path to the clone of the ighrealignment repo)
    *) usage ;;
  esac
done

# Clip path and suffixes from filename.
prefix=$(basename $bam)
prefix=${prefix/.sorted/}
prefix=${prefix/.deduped/}     # That's what the previous step of pipeline uses as suffix
prefix=${prefix/.bam/}
prefix=${prefix/.aln/}
#############

## Module Loading ##
module load bowtie/1.2.2_p1-foss-2017a
module load samtools/1.9-foss-2017a
module load python/2.7.13-foss-2017a
module load pysam/0.14.1-foss-2017a-python-2.7.13
module load pybedtools/0.8.0-foss-2017a-python-2.7.13
module load intervaltree/3.0.2-foss-2017a-python-2.7.13
#####################

# Create destination.
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
  mkdir -p $outDir
fi

echo "${prefix}: Sorting and indexing BAM"
samtools sort -@ 3 -m 9800M -o ${outDir}/${prefix}.sorted.bam $bam
samtools index -@ 3 ${outDir}/${prefix}.sorted.bam
bam="${outDir}/${prefix}.sorted.bam"

echo " "
echo "${prefix}: Extracting reads mapped to $CELLLINE"
samtools view -h $bam $CELLLINE > ${outDir}/${prefix}.sam

echo " "
echo "${prefix}: Categorise reads mapped to $CELLLINE"
# Look for multimapping to the natural chromosomal coordinates.
python2.7 ${SCRIPT_PATH}/DiscardNonIgHReads.py -c ${outDir}/${prefix}.sam -a $bam -g $CELLLINE
# Exclude unmapped reads (does that extra step save time?)
samtools view -h -F 4 $bam | awk -v chr="$CELLLINE" '{if($3 != chr){print $0}}' | samtools view -Sb - > ${outDir}/${prefix}.temp.bam
# Gather all qualifying reads.
samtools merge -@ 3 -f ${outDir}/${prefix}.temp2.bam ${outDir}/${prefix}.temp.bam ${outDir}/${prefix}_ighfiltered.bam ${outDir}/${prefix}_exonfiltered.bam ${outDir}/${prefix}_unique.bam
# Tidy
samtools sort -@ 3 -m 9800M -o ${outDir}/${prefix}_AllMapped.sorted.bam ${outDir}/${prefix}.temp2.bam
samtools index -@ 3 ${outDir}/${prefix}_AllMapped.sorted.bam
rm ${outDir}/${prefix}.temp.bam*
rm ${outDir}/${prefix}.temp2.bam*

echo ""
echo "${prefix}: IgH filtering complete"
rm ${outDir}/${prefix}.sorted.bam*

exit $?
#####################
