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
# Last reviewed on 05/aug/2019 by kimon.froussios@imp.ac.at
####################


set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -b BAM -o OUT_DIR -t PROtype -s CHRSIZES"
    exit 1
}
# Parse options.
while getopts 'b:o:t:s:' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;		     # Bamfile
    o) outDir="${OPTARG}" ;;		# output directory
    t) protype="${OPTARG}" ;;			# procap or proseq
    s) chrinfo="${OPTARG}" ;;			# chromosome sizes file
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


# Create destination.
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
  mkdir -p ${outDir}
fi

if [ "$PROprot" = "procap" ];then
	echo "${prefix}: Analyzing 5 prime end for PRO-cap analysis"
	NPRIME="5"
elif [ "$PROprot" = "proseq" ];then
	echo "${prefix}: Analyzing 3 prime end for PRO-seq analysis"
	NPRIME="3"
else
	exit
fi
#############

## Module Loading ##
module load samtools/1.9-foss-2017a
module load bedtools/2.27.1-foss-2017a
module load ucsc-kent-utils/356-foss-2017a
####################

#####################
mkdir -p ${outDir}
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
fi

echo ""
echo "${prefix}: Sorting bam by coordinates"
samtools sort -@3 -m 1000M -o ${outDir}/${prefix}.sorted.bam $bam

echo ""
echo "${prefix}: Getting mapped read counts."
AllMappedReads=$(samtools view -F 4 -c ${outDir}/${prefix}.sorted.bam)
#AllMappedReads=$(samtools view ${outDir}/${prefix}.sorted.bam | cut -f 1 | sort | uniq | wc -l)
echo ""
echo "${prefix}: Number of mappable reads: ${AllMappedReads}"
rpmFactor=$(awk -v mr="${AllMappedReads}" 'BEGIN {printf 1000000/mr}')
echo ""
echo "${prefix}: rpmFactor"
echo $rpmFactor

echo ""
echo "${prefix}: Assembling bedGraph from BAM."
bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -${NPRIME} -strand + -scale ${rpmFactor} > ${outDir}/${prefix}.${NPRIME}prime_rs_AllMappedReads_plus.unsorted.bedGraph
bedtools genomecov -ibam ${outDir}/${prefix}.sorted.bam -bg -${NPRIME} -strand - -scale ${rpmFactor} > ${outDir}/${prefix}.temp.bedGraph
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,-1*$4}' ${outDir}/${prefix}.temp.bedGraph > ${outDir}/${prefix}.${NPRIME}prime_rs_AllMappedReads_minus.unsorted.bedGraph
rm ${outDir}/${prefix}.temp.bedGraph

bedSort ${outDir}/${prefix}.${NPRIME}prime_rs_AllMappedReads_plus.unsorted.bedGraph ${outDir}/${prefix}_${NPRIME}prime_AllMappedReads_plus.bedGraph
bedSort ${outDir}/${prefix}.${NPRIME}prime_rs_AllMappedReads_minus.unsorted.bedGraph ${outDir}/${prefix}_${NPRIME}prime_AllMappedReads_minus.bedGraph
rm ${outDir}/${prefix}.${NPRIME}prime_rs_AllMappedReads_plus.unsorted.bedGraph
rm ${outDir}/${prefix}.${NPRIME}prime_rs_AllMappedReads_minus.unsorted.bedGraph

echo ""
echo "${prefix}: Converting bedGraph to BigWig"
bedGraphToBigWig ${outDir}/${prefix}_${NPRIME}prime_AllMappedReads_plus.bedGraph $chrinfo ${outDir}/${prefix}_${NPRIME}prime_AllMappedReads_plusstrand.bw
bedGraphToBigWig ${outDir}/${prefix}_${NPRIME}prime_AllMappedReads_minus.bedGraph $chrinfo ${outDir}/${prefix}_${NPRIME}prime_AllMappedReads_minusstrand.bw

echo " "
echo "${prefix}: BigWig files created"
rm ${outDir}/${prefix}.sorted.bam*

exit $?
#############
