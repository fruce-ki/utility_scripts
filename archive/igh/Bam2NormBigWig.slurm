#!/bin/bash
#
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=7000M
#SBATCH --cpus-per-task=4

#############
## Read Me ##
#############
# Goal of script: Convert Bam file to a Normalized bigwig format
# Pipeline based on the one proposed by Danko Lab: https://github.com/Danko-Lab/utils/blob/master/proseq/proseqMapper.bsh
#
# Last reviewed on 06/jul/2020 by kimon.froussios@imp.ac.at
####################


set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -b BAM -o OUT_DIR -t PROtype -s CHRSIZES -H METRICS -x SCRIPT_PATH [-L LIBSIZE | -F FACTOR]"
    exit 1
}
# Parse options.
scaleFactor=1
while getopts 'b:o:t:s:H:x:L:F:' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;		     # Bamfile
    o) outDir="${OPTARG}" ;;		# output directory
    t) protype="${OPTARG}" ;;		# procap or proseq or quantrev
    s) chrinfo="${OPTARG}" ;;		# chromosome sizes file
    H) metrics="${OPTARG}" ;;		# read counts
    x) scriptdir="${OPTARG}" ;;		# path to shiftBedgraphCoordinates.R (ie. the local clone of the IgH code repo)
    F) scaleFactor="${OPTARG}" ;;	# Scale values by this factor
    L) AllMappedReads="${OPTARG}" ;;	# Scale to RPM using this whole-library size for factor calculation
    *) usage ;;
  esac
done

if [ -z "$bam" ] || [ -z "$outDir" ] || [ -z "$protype" ]; then
	echo "${prefix}: Insufficient parameters."
	usage
	exit 1
fi

PROprot=$(echo "${protype}" | awk '{print tolower($0)}')
GENOME=$(echo "${GENOME}" | awk '{print tolower($0)}')
#prmtr="m60_mm"
####################

prefix=$(basename $bam)
prefix=${prefix/.sorted/}
prefix=${prefix/.bam/}
prefix=${prefix/.aln/}
prefix=${prefix/.deduped/}


if [ "$PROprot" = "procap" ];then
	NPRIME="5"
elif [ "$PROprot" = "proseq" ];then
	NPRIME="3"
elif [ "$PROprot" = "quantrev" ];then
    NPRIME="3"
else
	exit
fi
#############

## Module Loading ##
# module load samtools/1.9-foss-2018b
# module load bedtools/2.27.1-foss-2018b
# module load kent_tools/20190507-linux.x86_64
####################

#####################

if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
    mkdir -p ${outDir}
fi

echo ""
echo "${prefix}: Sorting bam by coordinates"
samtools sort -@3 -m 1500M -o ${outDir}/${prefix}.sorted.bam $bam
bam=${outDir}/${prefix}.sorted.bam

# echo ""
# echo "${prefix}: Getting mapped read counts."
# # AllMappedReads=$(samtools view ${outDir}/${prefix}.sorted.bam | cut -f 1 | sort | uniq | wc -l)
# # AllMappedReads=$(samtools view -F 4 -c ${bam})
# trackReads=$(srun ${scriptdir}/CountMappedReadsBAM.py -b $bam | cut -f 3)
# echo "${prefix} STATS: $trackReads reads for the track"
# printf '%s\t%d\n' $(perl -e '$ARGV[0]=~/(unique|nonunique|locus)/; print $1' $prefix) $trackReads >> $metrics

# Calculate RPM factor
if [ ! -z "$AllMappedReads" ]; then
    echo ""
    scaleFactor=$(awk -v mr="${AllMappedReads}" 'BEGIN {printf 1000000/mr}')
    echo "${prefix}: RPM scaling factor = $scaleFactor"
elif (( $(echo "$scaleFactor > 1" | bc -l) )); then
    echo "${prefix}: custom scaling factor = $scaleFactor"
fi

echo ""
echo "${prefix}: Assembling bedGraph from BAM."
if [ "$PROprot" == 'procap' ]; then
    echo "${prefix}: PROcap stranding"
    # keeping strand designations (PROcap reads are the same strand as the mRNA)
    # 5' of reads == 5' of mRNA
    bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -5 -strand + -scale ${scaleFactor} > ${outDir}/${prefix}.${NPRIME}prime+.unsorted.bedGraph
    bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -5 -strand - -scale ${scaleFactor} > ${outDir}/${prefix}.temp.bedGraph
elif [ "$PROprot" == 'proseq' ]; then
    echo "${prefix}: PROseq stranding"
    # flipping strand designations (PROseq reads are the opposite strand from the mRNA). This way PROseq and PROcap tracks both match the mRNA strand.
    # 5' of the reads == 3' of mRNA
    bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -5 -strand - -scale ${scaleFactor} > ${outDir}/${prefix}.${NPRIME}prime+.unsorted.bedGraph
    bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -5 -strand + -scale ${scaleFactor} > ${outDir}/${prefix}.temp.bedGraph
    # Shift coordinates by 1 (the 2nd base of the read is the last "real" base of the mRNA, before intervening with the termination NTPs)
    # Also remember that the strand labels are artificially flipped. So the "+" are actually mapped to - and need to be shifted left, and the "-" are actually mapped to + and need to be shifted right.
    ${scriptdir}/shiftBedgraphCoordinates.R ${outDir}/${prefix}.${NPRIME}prime+.unsorted.bedGraph '-1' ${outDir}/${prefix}.${NPRIME}prime+.unsorted.bedGraph
    ${scriptdir}/shiftBedgraphCoordinates.R ${outDir}/${prefix}.temp.bedGraph '1' ${outDir}/${prefix}.temp.bedGraph
elif [ "$PROprot" == 'quantrev' ]; then
    echo "${prefix}: quantseq-REV stranding"
    # flipping strand designations (Quantseq-REV reads are the opposite strand from the mRNA).
    # 5' of the reads == 3' of mRNA
    bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -5 -strand - -scale ${scaleFactor} > ${outDir}/${prefix}.${NPRIME}prime+.unsorted.bedGraph
    bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -5 -strand + -scale ${scaleFactor} > ${outDir}/${prefix}.temp.bedGraph
    # No need to shift coordinates. Supposedly the first base of quantseq-REV reads is the last real base of the mRNA before the polyA tail.
fi
# Multiply - strand scores by -1, for nicer visualisation in genome browser.
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,-1*$4}' ${outDir}/${prefix}.temp.bedGraph > ${outDir}/${prefix}.${NPRIME}prime-.unsorted.bedGraph
rm ${outDir}/${prefix}.temp.bedGraph

bedSort ${outDir}/${prefix}.${NPRIME}prime+.unsorted.bedGraph ${outDir}/${prefix}_${NPRIME}prime+.bedGraph
bedSort ${outDir}/${prefix}.${NPRIME}prime-.unsorted.bedGraph ${outDir}/${prefix}_${NPRIME}prime-.bedGraph
rm ${outDir}/${prefix}.${NPRIME}prime+.unsorted.bedGraph
rm ${outDir}/${prefix}.${NPRIME}prime-.unsorted.bedGraph

echo ""
echo "${prefix}: Converting bedGraph to BigWig"
bedGraphToBigWig ${outDir}/${prefix}_${NPRIME}prime+.bedGraph $chrinfo ${outDir}/${prefix}_${NPRIME}prime+.bw || true
bedGraphToBigWig ${outDir}/${prefix}_${NPRIME}prime-.bedGraph $chrinfo ${outDir}/${prefix}_${NPRIME}prime-.bw || true
# Conversion fails if bedgraph empty, which is fine, but I still want the sorted bam to be removed.

echo " "
echo "${prefix}: Tracks created"
rm ${outDir}/${prefix}.sorted.bam*

exit $?
#############
