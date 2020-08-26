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
# Last reviewed: 07/jan/2020	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -b BAM_FILE -o OUTDIR -c CELLLINE -H METRICS [-x SCRIPTDIR]"
    exit 1
}
# Defaults.
SCRIPT_PATH='/groups/pavri/Kimon/ighrealignment'
# Parse options.
while getopts 'b:o:c:l:x:H:' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;           # BAM file.
    o) outDir="${OPTARG}" ;;        # Output directory.
    c) CELLLINE="${OPTARG}" ;;      # Full reference name of locus chromosome.
    x) SCRIPT_PATH="${OPTARG}" ;;   # Directory with the python scripts (ie. path to the clone of the ighrealignment repo)
    H) metrics="${OPTARG}" ;;		# read counts
    *) usage ;;
  esac
done

# Clip path and possible suffixes from filename.
prefix=$(basename $bam)
prefix=${prefix/.sorted/}
prefix=${prefix/.deduped/}
prefix=${prefix/.aln/}
prefix=${prefix/_subject/}
prefix=${prefix/.bam/}
prefix=${prefix/.aln/}
#############

## Module Loading ##
# module load samtools/1.9-foss-2018b
# module load python/3.6.6-foss-2018b
# module load pysam/0.15.1-foss-2018b-python-3.6.6
# module load pybedtools/0.8.0-foss-2018b-python-3.6.6
# module load intervaltree/3.0.2-foss-2018b-python-3.6.6
#####################

# Create destination.
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
  mkdir -p $outDir
fi

# cnt=$(samtools view -c $bam)
cnt=$(${SCRIPT_PATH}/CountMappedReadsBAM.py -b $bam | cut -f 3)   # Bowtie does not mark secondary alignments, so samtools count is inflated by multimappers
echo "${prefix} STATS: $cnt unfiltered reads"
printf '%s\t%d\n' "unfiltered" $cnt >> $metrics

echo ""
echo "${prefix}: Categorize reads mapped to ${CELLLINE}"
${SCRIPT_PATH}/DiscardNonIgHReads.py -a $bam -g $CELLLINE -o ${outDir}

unique="${outDir}/${prefix}_unique.bam"
nonunique="${outDir}/${prefix}_nonunique.bam"
target="${outDir}/${prefix}_locus.bam"

# cnt=$(samtools view -c $target)
cnt=$(${SCRIPT_PATH}/CountMappedReadsBAM.py -b $target | cut -f 3)
echo "${prefix} STATS: $cnt reads mapping to locus"
printf '%s\t%d\n' "_locus" $cnt >> $metrics

# cnt=$(samtools view -c $unique)
cnt=$(${SCRIPT_PATH}/CountMappedReadsBAM.py -b $unique | cut -f 3)
echo "${prefix} STATS: $cnt reads mapping uniquely to locus"
printf '%s\t%d\n' "__unique" $cnt >> $metrics

# cnt=$(samtools view -c $nonunique)
cnt=$(${SCRIPT_PATH}/CountMappedReadsBAM.py -b $nonunique | cut -f 3)
echo "${prefix} STATS: $cnt reads multi-mapping to locus"
printf '%s\t%d\n' "__nonunique" $cnt >> $metrics

echo ""
echo "${prefix}: IgH filtering complete"

exit $?
#####################
