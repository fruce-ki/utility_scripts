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
# Last reviewed: 07/jan/2020	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      Bam2Fq -i BAM_FILE -o OUT_DIR -H METRICS -a 3_ADAPTOR -b 5_ADAPTOR [-l MIN_LENGTH] [-u umiLen] [-q MIN_Q]"
    exit 1
}
# Defaults

MinLen=12
MinQ=30
umiLen=0
# Parse options.
while getopts 'i:o:a:l:q:u:H:b:' flag; do
  case "${flag}" in
    i) bam="${OPTARG}" ;;			# BAM file
    o) outDir="${OPTARG}" ;;		# output directory.
    a) adapt3="${OPTARG}" ;;		# 3' adaptor, including the -a flag for cutadapt, ie.  "-a ATGC"
    b) adapt5="${OPTARG}" ;;		# 5' adaptor, including the -g flag for cutadapt, ie.  "-g ATGC"
    l) MinLen="${OPTARG}" ;;		# minimum acceptable read length (12)
    q) MinQ="${OPTARG}" ;;		    # minimum acceptable base quality at 3' (30)
    H) metrics="${OPTARG}" ;;		# read counts
    u) umiLen="${OPTARG}" ;;		# length of UMIs (0) (will be added to MinLen for the actual minimum trimmed length)
    *) usage ;;
  esac
done

if [ ! -z "$adapt3" ]; then
    adapt3="-a ${adapt3}"
fi
if [ ! -z "$adapt5" ]; then
    adapt5="-g ${adapt5}"
fi

# Clip path and suffixes from filename.
prefix=$(basename $bam)
prefix=${prefix/.bam/''}
prefix=${prefix/.sorted/''}
prefix=$(perl -e 'if($ARGV[0]=~/(.+)_\d{8}\w?_\d{8}$/){print $1}else{print $ARGV[0]}' $prefix)   # Get rid of the date, to make names shorter and keep them manageable when more things are appended to them.
####################

## Module Loading ##
# module load cutadapt/1.18-foss-2018b-python-3.6.6
# module load samtools/1.9-foss-2018b
#####################

# Create destination.
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
  mkdir -p ${outDir} ${outDir}_qc ${outDir}_qc_trimmed
fi

# # If BAM not sorted, sort it.
# cleanup=0
# samtools view -H $bam | head -n1 > ${outDir}/${prefix}.temp_header.txt
# if grep -E -q "SO:coordinate" ${outDir}/${prefix}.temp_header.txt; then
# 	echo "${prefix}: sorted bam as input"
# else
# 	echo "${prefix}: sorting input bam"
# 	samtools sort -@ 3 -o ${outDir}/${prefix}.sorted.bam $bam
# 	bam="${outDir}/${prefix}.sorted.bam"
#   cleanup=1
# fi
# rm ${outDir}/${prefix}.temp_header.txt

cnt=$(samtools view -c ${bam})
echo "${prefix} STATS: $cnt unprocessed reads"
printf '%s\t%d\n' "unprocessed" $cnt >> $metrics

echo "${prefix}: Quality report"
fastqc -t 4 -o ${outDir}_qc -f bam $bam


echo "${prefix}: Extracting reads"
samtools bam2fq -@ 3 -0 ${outDir}/${prefix}.MyFq.fq ${bam}
fq="${outDir}/${prefix}.MyFq.fq"


echo "${prefix}: Trimming 5'/3' adaptors and 3' qualities <${MinQ} & filtering reads <${MinLen}"
MinLen=$((MinLen + umiLen))
# fastx_clipper -d200 -Q33 -l${MinLen} -i ${outDir}/${prefix}.MyFq.fq -o ${outDir}/${prefix}.fq
cutadapt -j 4 -q $MinQ -m $MinLen $adapt3 $adapt5 --report full -o ${outDir}/${prefix}.trimmed.fq $fq > ${outDir}/${prefix}.cutadapt.log
rm $fq
fq=${outDir}/${prefix}.trimmed.fq

cnt=$(echo "$(wc -l <$fq)/4" | bc)
echo "${prefix} STATS: $cnt trimmed/filtered reads"
printf '%s\t%d\n' "trimmed" $cnt >> $metrics

echo "${prefix}: Quality report"
fastqc -t 4 -o ${outDir}_qc_trimmed -f fastq $fq


echo "${prefix}: Compressing FASTQ file"
gzip -f $fq


echo "${prefix}: Finished read pre-processing."
exit $?
####################
