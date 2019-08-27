#!/bin/bash

#############
## Read Me ##
#############
# A meta-script running the main IgH realignement steps.
#
# Last reviewed: 05/aug/2019	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      igh_workflow.slurm -i BAM_FILE -o OUT_DIR -L LOCUS -t PROtype [-m MISMATCHES] [-x SCRIPT_DIR] [-n UMI_LENGTH] [-u] [-k]"
    exit 1
}
# Defaults.
SCRIPT_PATH='/groups/pavri/Kimon/ighrealignment'
mm_tol=3
spikedin=0
clipped=0
# Parse options.
while getopts 'i:o:L:t:m:n:x:uk' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;				# BAM file for one sample.
    o) outDir="${OPTARG}" ;;		# output root directory. fastq, alignment, filtering and tracks subdirectories will be created in there)
    L) LOCUS="${OPTARG}" ;;     # Predefined options for index.
    t) protype="${OPTARG}" ;;   # procap or proseq
    m) mm_tol="${OPTARG}" ;;    # number of mismatches allowed
    n) umi="${OPTARG}" ;;       # UMI length. Use 0 to skip deduplication steps.
    x) SCRIPT_PATH="${OPTARG}" ;;	 # Directory with the python scripts (ie. path to the clone of the ighrealignment repo)
    k) spikedin=1 ;;              # Use index that includes the dm_r6 genome.
    u) clipped=1 ;;              # UMIs are already clipped and added to end of read title.
    *) usage ;;
  esac
done
####################


#############
# Clip suffixes from sample.
prefix=$(basename $bam)
prefix=${prefix/.bam/''}
prefix=${prefix/.sorted/''}
prefix=$(perl -e 'if($ARGV[0]=~/(.+)_\d{8}\w?_\d{8}$/){print $1}else{print $ARGV[0]}' $prefix)   # Get rid of the date, to make names shorter and keep them manageable when more things are appended to them.

echo "Output prefix: ${prefix}"

# Save command-line typing and mistakes by pre-defining the index.
## !!! When adding entries here, remember to also: !!!
##      - update the IntervalTree definitions in DiscardNonIgHReads.py
##      - update the chromosome sizes file or defining a new one
if [ "$LOCUS" = "B18" ]; then
  INDEXDIR="/users/kimon.froussios/ursi/PRO/aux/bowtie"
  CELLLINE="B1-8hi_mm10_190516_181014-most-recent"
  CHRINFO="/users/kimon.froussios/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt" # file containing the chromosome sizes of the genome used
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDRc" ]; then
  INDEXDIR="/users/kimon.froussios/ursi/PRO/aux/bowtie"
	CELLLINE="B18hi_HDRc_BfaI_ctrl_190513_US"
  CHRINFO="/users/kimon.froussios/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR1" ]; then
  INDEXDIR="/users/kimon.froussios/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR1_dTATA_190513_US"
  CHRINFO="/users/kimon.froussios/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR2" ]; then
  INDEXDIR="/users/kimon.froussios/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR2_dTATA_+1T_190513_US"
  CHRINFO="/users/kimon.froussios/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR3" ]; then
  INDEXDIR="/users/kimon.froussios/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR3_TATTAC_190513_US"
  CHRINFO="/users/kimon.froussios/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR4" ]; then
  INDEXDIR="/users/kimon.froussios/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR4_TATTAC_+1T_190513_US"
  CHRINFO="/users/kimon.froussios/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
else
  echo "${prefix}: Invalid predefined LOCUS ${LOCUS}"
  exit 1
fi

# Did the samples include spike-in of Drosophila melanogaster nuclei?
if [ $spikedin -eq 1 ]; then
  INDEX="${INDEXDIR}/${GENOME}_plus_dmr6_plus_${CELLLINE}"
else
  INDEX="${INDEXDIR}/${GENOME}_plus_${CELLLINE}"
fi

echo ""
echo "${prefix}: BAM2FASTQ"
srun --qos short --cpus-per-task=4 --ntasks=1 --nodes=1 -J $prefix ${SCRIPT_PATH}/Bam2Fq.slurm -i $bam -o ${outDir}/fastq

echo ""
echo "${prefix}: ALIGN_&_DEDUPLICATE"
srun --qos short --mem=40000M --cpus-per-task=4 --ntasks=1 --nodes=1 -J $prefix ${SCRIPT_PATH}/Align_and_Deduplicate.slurm  -f ${outDir}/fastq/${prefix}.fq.gz -o ${outDir}/alignment -i $INDEX -r $CHRINFO -l $LOCUS -m $mm_tol -n $umi -u $clipped

echo ""
echo "${prefix}: FILTER_IgH"
if [ -f "${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.deduped.sorted.bam" ]; then
  srun --qos short --mem=40000M --cpus-per-task=4 --ntasks=1 --nodes=1 -J $prefix ${SCRIPT_PATH}/Filter_IgH.slurm -b ${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.deduped.sorted.bam -o ${outDir}/filtering -c $CELLLINE -x $SCRIPT_PATH
else
  srun --qos short --mem=40000M --cpus-per-task=4 --ntasks=1 --nodes=1 -J $prefix ${SCRIPT_PATH}/Filter_IgH.slurm -b ${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.aln.sorted.bam -o ${outDir}/filtering -c $CELLLINE -x $SCRIPT_PATH
fi

echo ""
echo "${prefix}: TRACKS"
srun --qos=short --cpus-per-task=4 --mem=5000M --ntasks=1 --nodes=1 -J $prefix ${SCRIPT_PATH}/Bam2NormBigWig.slurm -b ${outDir}/filtering/${prefix}_${LOCUS}_mm${mm_tol}_AllMapped.sorted.bam -o ${outDir}/tracks -t $protype -s $CHRINFO

echo ""
echo "${prefix}: WORKFLOW FINISHED"
#############
