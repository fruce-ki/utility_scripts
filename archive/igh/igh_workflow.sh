#!/bin/bash

#############
## Read Me ##
#############
# A meta-script running the main IgH realignement steps.
#
# Last reviewed: 07/jan/2020	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      igh_workflow.slurm {-1 -2 -3 -4 -5} -i BAM_FILE -o OUT_DIR -L LOCUS -t PROtype [-m MISMATCHES] [-x SCRIPT_DIR] [-n UMI_LENGTH] [-u] [-U] [-k]"
    exit 1
}
# Defaults.
SCRIPT_PATH='/groups/pavri/Kimon/ighrealignment'
mm_tol=3
spikedin=0
splitnorm=0
umi=0
clipped=0
nodedup=0
bu2fq=0
fq2ba=0
ba2fl=0
fl2wig=0
multis=194                              # (chosen because there are 183+11 IgH V-segments in mm10.6. Adjust as appropriate)
#PROseq
proseq5="GTTCAGAGTTCTACAGTCCGACGATC"        # 3' mRNA adaptor rev'comp'ed for read
proseq3="NNNNTGGAATTCTCGGGTGCCAAGG"         # 5' mRNA adaptor rev'comp'ed for read
#PROcap
procap3="NNNNAGATCGGAAGAGCACACGTCT"            # 3' mRNA and read adaptor
procap5="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"    # 5' mRNA and read adaptor
#quantrev
quantrev5="TTTTTTTTTTTTTTTTTTT"                # polyA revcomp

# Parse options.
while getopts 'i:o:L:t:m:n:x:M:uUk12345' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;				# BAM file for one sample.
    o) outDir="${OPTARG}" ;;		# output root directory. fastq, alignment, filtering and tracks subdirectories will be created in there)
    L) LOCUS="${OPTARG}" ;;     # Predefined options for index.
    t) protype="${OPTARG}" ;;   # procap, proseq or quantrev
    m) mm_tol="${OPTARG}" ;;    # number of mismatches allowed
    n) umi="${OPTARG}" ;;       # UMI length (0)
    x) SCRIPT_PATH="${OPTARG}" ;;	 # Directory with the python scripts (ie. path to the clone of the ighrealignment repo)
    M) multis="${OPTARG}" ;;		# Don't report multimappers with more mappings than this many valid alignments (194).
    k) spikedin=1 ;;              # Use index that includes the dm_r6 genome.
    u) clipped=1 ;;              # UMIs are already clipped and added to end of read title.
    U) nodedup=1 ;;                # 1: Clip the UMI, but don't deduplicate.
    1) bu2fq=1 ;;                # Do Bam2Fq
    2) fq2ba=1 ;;                # Do align and deduplicate
    3) splitnorm=1 ;;            # Split reads by species
    4) ba2fl=1 ;;                # Do re-assignment
    5) fl2wig=1 ;;               # Do tracks
    *) usage ;;
  esac
done

wait_for_jobs(){
  echo "waiting"
  sleep 60  # seconds, give time to the scheduler to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep kimon.fr | grep -c $1) -eq 0 ]; then
      break
    else
      echo "sleep another" $((sleeptime / 60)) "minutes..."
      sleep $sleeptime
    fi
    n=$((n + 1))
    if [ $n -eq 5 ]; then
      sleeptime=300  # if still running after 10 mins, ask every 5 mins
    fi
    if [ $n -eq 10 ]; then
      sleeptime=600  # if still running after 30 mins, ask every 10 mins
    fi
  done
}
####################


#############
# Clip suffixes from sample.
prefix=$(basename $bam)
prefix=${prefix/.bam/''}
prefix=${prefix/.sorted/''}
prefix=$(perl -e 'if($ARGV[0]=~/(.+)_\d{8}\w?_\d{8}$/){print $1}else{print $ARGV[0]}' $prefix)   # Get rid of the date, to make names shorter and keep them manageable when more things are appended to them.


protype=$(echo "${protype}" | awk '{print tolower($0)}')


echo "Output prefix: ${prefix}"

# Save command-line typing and mistakes by pre-defining the index.
## !!! When ADDING ENTRIES here, REMEMBER to also:
##      - update the IntervalTree definitions in DiscardNonIgHReads.py
##      - update the chromosome sizes file or define a new one
if [ "$LOCUS" = "B18" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/PRO/aux/bowtie"
  CELLLINE="B1-8hi_mm10_190516_181014-most-recent"
  CHRINFO="/groups/pavri/Kimon/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt" # file containing the chromosome sizes of the genome used
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDRc" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/PRO/aux/bowtie"
	CELLLINE="B18hi_HDRc_BfaI_ctrl_190513_US"
  CHRINFO="/groups/pavri/Kimon/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR1" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR1_dTATA_190513_US"
  CHRINFO="/groups/pavri/Kimon/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR2" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR2_dTATA_+1T_190513_US"
  CHRINFO="/groups/pavri/Kimon/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR3" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR3_TATTAC_190513_US"
  CHRINFO="/groups/pavri/Kimon/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "B18_HDR4" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/PRO/aux/bowtie"
  CELLLINE="B18hi_HDR4_TATTAC_+1T_190513_US"
  CHRINFO="/groups/pavri/Kimon/ursi/PRO/aux/mm10-dmr6-B18_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "CH12" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/CH12_PRO/aux/bowtie1"
  CELLLINE="CH12_VDJ_200122"
  CHRINFO="/groups/pavri/Kimon/ursi/CH12_PRO/aux/mm10-dmr6-CH12_chrSizes.txt"
  GENOME="mm10"              # NCBI mm10.6
#  CHROMOSOME="NC_000078.6"  # NCBI chr12
elif [ "$LOCUS" = "Ramos_IGH" ]; then
  INDEXDIR="/groups/pavri/Kimon/ursi/Ramos_PRO/aux/Hg38_plus_RamosIgH"
  CELLLINE="Ramos_IgH"
  CHRINFO="/groups/pavri/Kimon/ursi/Ramos_PRO/aux/Hg38_plus_RamosIgH/Hg38_plus_RamosIgH_chrsizes.txt"
  GENOME="Hg38"
else
  echo "${prefix}: Invalid predefined LOCUS ${LOCUS}"
  exit 1
fi

# Did the samples include spike-in of Drosophila melanogaster nuclei?
if [ "$spikedin" -eq 1 ]; then
  INDEX="${INDEXDIR}/${GENOME}_plus_dmr6_plus_${CELLLINE}"
else
  INDEX="${INDEXDIR}/${GENOME}_plus_${CELLLINE}"
fi

if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}
fi


# Metrics
metrics="${outDir}/${prefix}.counts.txt"
printf '%s\t%s\n' "Stage" "Reads" >> $metrics

if [ "$bu2fq" -eq 1 ]; then
    echo ""
    echo "${prefix}: BAM2FASTQ"
    if [ "$protype" == 'procap' ]; then
        adapt3="$procap3"
        adapt5="$procap5"
    elif [ "$protype" == 'proseq' ]; then
        adapt3="$proseq3"
        adapt5="$proseq5"
    elif [ "$protype" == 'quantrev' ]; then
        adapt3=""
        adapt5="$quantrev5"
    fi

    if [ "$protype" == 'quantrev' ]; then
        # -o /dev/null -e /dev/null
        sbatch --qos short --cpus-per-task=4 --ntasks=1 --nodes=1 -J IGHb2fq --wrap "${SCRIPT_PATH}/Bam2Fq.slurm -i $bam -o ${outDir}/fastq -u $umi -H $metrics"
    else
        # -o /dev/null -e /dev/null
        sbatch --qos short --cpus-per-task=4 --ntasks=1 --nodes=1 -J IGHb2fq --wrap "${SCRIPT_PATH}/Bam2Fq.slurm -i $bam -o ${outDir}/fastq -u $umi -a $adapt3 -b $adapt5 -H $metrics"
    fi
    wait_for_jobs IGHb2fq
fi

if [ "$fq2ba" -eq 1 ]; then
    echo ""
    echo "${prefix}: ALIGN_&_DEDUPLICATE"
    # -o /dev/null -e /dev/null
    sbatch --qos medium --mem=40000M --cpus-per-task=4 --ntasks=1 --nodes=1 -J IGHaln --wrap "${SCRIPT_PATH}/Align_and_Deduplicate.slurm  -f ${outDir}/fastq/${prefix}.trimmed.fq.gz -o ${outDir}/alignment -i $INDEX -r $CHRINFO -l $LOCUS -m $mm_tol -n $umi -u $clipped -U $nodedup -M $multis -H $metrics -x $SCRIPT_PATH"
    wait_for_jobs IGHaln
fi

if [ "$splitnorm" -eq 1 ] && [ "$spikedin" -eq 1 ]; then
    echo ""
    echo "${prefix}: SPIKE_INs"
    # Depending on whether deduplication was applied:
    if [ -f "${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam" ]; then
        # -o /dev/null -e /dev/null
        sbatch --qos short --mem=40000M --cpus-per-task=1 --ntasks=1 --nodes=1 -J IGHspk --wrap "${SCRIPT_PATH}/Filter_SpikeIn.slurm -b ${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam -o ${outDir}/spikein -x $SCRIPT_PATH -H $metrics"
    else
        # -o /dev/null -e /dev/null
        sbatch --qos short --mem=40000M --cpus-per-task=1 --ntasks=1 --nodes=1 -J IGHspk --wrap "${SCRIPT_PATH}/Filter_SpikeIn.slurm -b ${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.aln.bam -o ${outDir}/spikein -x $SCRIPT_PATH -H $metrics"
    fi
    wait_for_jobs IGHspk
fi

if [ "$ba2fl" -eq 1 ]; then
    echo ""
    echo "${prefix}: FILTER_IgH"
    # Depending on whether deduplication and spike-in separation were applied:
    if [ -f "${outDir}/spikein/${prefix}_${LOCUS}_mm${mm_tol}_subject.bam" ]; then
        # -o /dev/null -e /dev/null
        sbatch --qos short --mem=40000M --cpus-per-task=4 --ntasks=1 --nodes=1 -J IGHfltr --wrap "${SCRIPT_PATH}/Filter_IgH.slurm -b ${outDir}/spikein/${prefix}_${LOCUS}_mm${mm_tol}_subject.bam -o ${outDir}/filtering -c $CELLLINE -x $SCRIPT_PATH -H $metrics"
    elif [ -f "${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam" ]; then
        # -o /dev/null -e /dev/null
        sbatch --qos short --mem=40000M --cpus-per-task=4 --ntasks=1 --nodes=1 -J IGHfltr --wrap "${SCRIPT_PATH}/Filter_IgH.slurm -b ${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam -o ${outDir}/filtering -c $CELLLINE -x $SCRIPT_PATH -H $metrics"
    else
        # -o /dev/null -e /dev/null
        sbatch --qos short --mem=40000M --cpus-per-task=4 --ntasks=1 --nodes=1 -J IGHfltr --wrap "${SCRIPT_PATH}/Filter_IgH.slurm -b ${outDir}/alignment/${prefix}_${LOCUS}_mm${mm_tol}.aln.bam -o ${outDir}/filtering -c $CELLLINE -x $SCRIPT_PATH -H $metrics"
    fi
    wait_for_jobs IGHfltr
fi

if [ "$fl2wig" -eq 1 ]; then
    echo ""
    echo "${prefix}: TRACKS"
    echo "${outDir}/filtering/${prefix}_${LOCUS}_mm${mm_tol}_locus.bam"
    # -o /dev/null -e /dev/null
    sbatch --qos=short --cpus-per-task=4 --mem=5000M --ntasks=1 --nodes=1 -J IGHbed --wrap "${SCRIPT_PATH}/Bam2NormBigWig.slurm -b ${outDir}/filtering/${prefix}_${LOCUS}_mm${mm_tol}_locus.bam -o ${outDir}/tracks -t $protype -s $CHRINFO -H $metrics -x $SCRIPT_PATH"
    # -o /dev/null -e /dev/null
    sbatch --qos=short --cpus-per-task=4 --mem=5000M --ntasks=1 --nodes=1 -J IGHbed --wrap "${SCRIPT_PATH}/Bam2NormBigWig.slurm -b ${outDir}/filtering/${prefix}_${LOCUS}_mm${mm_tol}_nonunique.bam -o ${outDir}/tracks -t $protype -s $CHRINFO -H $metrics -x $SCRIPT_PATH"
    # -o /dev/null -e /dev/null
    sbatch --qos=short --cpus-per-task=4 --mem=5000M --ntasks=1 --nodes=1 -J IGHbed --wrap "${SCRIPT_PATH}/Bam2NormBigWig.slurm -b ${outDir}/filtering/${prefix}_${LOCUS}_mm${mm_tol}_unique.bam -o ${outDir}/tracks -t $protype -s $CHRINFO -H $metrics -x $SCRIPT_PATH"
    wait_for_jobs IGHbed
fi

# rm ./slurm*

echo ""
echo "${prefix}: WORKFLOW FINISHED"
#############
