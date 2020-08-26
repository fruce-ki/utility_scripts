#!/bin/bash
#
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40000

#############
## Read Me ##
#############
# It takes as input a (deduplicated) aligned BAM and the name of the dir where output files will be created.
#
# Created: 07/jan/2020	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -b BAM_FILE -o OUTDIR [-x SCRIPTDIR]"
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

# Clip path and suffixes from filename.
prefix=$(basename $bam)
prefix=${prefix/.sorted/}
prefix=${prefix/.deduped/}
prefix=${prefix/.bam/}
prefix=${prefix/.aln/}
#############

## Module Loading ##
# module load samtools/1.9-foss-2018b
# module load python/3.6.6-foss-2018b
# module load pysam/0.15.1-foss-2018b-python-3.6.6
#####################

# Create destination.
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
  mkdir -p $outDir
fi

# cnt=$(samtools view -c $bam)
${SCRIPT_PATH}/CountMappedReadsBAM.py -b $bam > ${outDir}/${prefix}.cnt   # Bowtie does not mark secondary alignments, so samtools count is inflated by multimappers
cnt=$(cut -f 3 ${outDir}/${prefix}.cnt)
cntm=$(cut -f 5 ${outDir}/${prefix}.cnt)
echo "${prefix} STATS: $cnt pre-split reads of which ${cntm} multimappers"
printf '%s\t%d\n' "mixed" "$cnt" >> $metrics
printf '%s\t%d\n' "_mixedmulti" "$cntm" >> $metrics

echo ""
echo "${prefix}: Split reads between spike-in and subject"
${SCRIPT_PATH}/SplitSpikein.py -a $bam -o ${outDir}


subject="${outDir}/${prefix}_subject.bam"
spikein="${outDir}/${prefix}_spikein.bam"
ambiguous="${outDir}/${prefix}_ambiguous.bam"

# cnt=$(samtools view -c $subject)
${SCRIPT_PATH}/CountMappedReadsBAM.py -b $subject > ${outDir}/${prefix}.cnt
cnt=$(cut -f 3 ${outDir}/${prefix}.cnt)
cntm=$(cut -f 5 ${outDir}/${prefix}.cnt)
echo "${prefix} STATS: $cnt reads mapping to subject of which ${cntm} multimappers"
printf '%s\t%d\n' "_subject" $cnt >> $metrics
printf '%s\t%d\n' "__subjectmutli" $cntm >> $metrics
rm ${outDir}/${prefix}.cnt

# cnt=$(samtools view -c $spikein)
${SCRIPT_PATH}/CountMappedReadsBAM.py -b $spikein > ${outDir}/${prefix}.cnt
cnt=$(cut -f 3 ${outDir}/${prefix}.cnt)
cntm=$(cut -f 5 ${outDir}/${prefix}.cnt)
echo "${prefix} STATS: $cnt reads mapping to spike-in of which ${cntm} multimappers"
printf '%s\t%d\n' "_spikein" $cnt >> $metrics
printf '%s\t%d\n' "__spikeinmutli" $cntm >> $metrics
rm ${outDir}/${prefix}.cnt

# cnt=$(samtools view -c $ambiguous)
cnt=$(${SCRIPT_PATH}/CountMappedReadsBAM.py -b $ambiguous | cut -f 3)
echo "${prefix} STATS: $cnt reads mapping to locus and IgH"
printf '%s\t%d\n' "_ambiguous" $cnt >> $metrics


echo ""
echo "${prefix}: Spike-in filtering complete"

exit $?
#####################
