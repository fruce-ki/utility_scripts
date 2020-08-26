#!/bin/bash
#
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G

#############
## Read Me ##
#############
# It takes as input a SINGLE-ENDED fastq file, the name of the dir where output files will be created,
# and a precomputed bowtie index.
#
# Last reviewed: 09/oct/2019	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -f FASTQ_FILE -o OUTDIR -i INDEX -r CHR_SIZES -l LOCUS -H METRICS -X SCRIPT_DIR [-m MISMATCHES] [-n UMI_LENGTH] [-u VAL] [-U VAL]"
    exit 1
}
# Defaults.
mm_tol=1
ntnum=10
mm_umi=0
clipped=0
nodedup=0
maxmultis=194   # (chosen because there are 183+11 IgH V-segments in mm10)
# Parse options.
while getopts 'f:o:i:l:r:m:n:M:H:u:U:x:' flag; do
  case "${flag}" in
    f) fqgz="${OPTARG}" ;;			# Compressed FASTQ file.
    o) outDir="${OPTARG}" ;;    # Output directory.
    i) INDEX="${OPTARG}" ;;			# Bowtie index prefix. (not bowtie2)
    r) ref_sizes="${OPTARGS}" ;;  # Table of reference sequence sizes
    l) LOCUS="${OPTARG}" ;;			# Short nickname for the index used, to be integrated into output file names.
    m) mm_tol="${OPTARG}" ;;		# Mismatch tolerance. Choose between 0 and 3 (1).
    M) maxmultis="${OPTARG}" ;;		# Ignore reads with more than this many valid alignments (194)
    n) ntnum="${OPTARG}" ;;			# Length of UMI at the start of reads (10). 0: No UMI.
    u) clipped="${OPTARG}" ;;       # 1: UMI already clipped and appended to end of read title, 0: UMI not clipped
    U) nodedup="${OPTARG}" ;;       # 1: Clip the UMI, but don't deduplicate. 0: clip and deduplicate
    H) metrics="${OPTARG}" ;;		# read counts
    x) XDIR="${OPTARG}" ;;
    *) usage ;;
  esac
done

if [ ${mm_tol} -eq 0 ] || [ ${mm_tol} -eq 1 ] || [ ${mm_tol} -eq 2 ] || [ ${mm_tol} -eq 3 ]; then
	echo ""
	echo "${prefix}: setting mismatch tolerance to ${mm_tol}"
else
	echo "${prefix}: Incorrect parameters! Mismatch tolerance has to be between 0 to 3"
	exit 1
fi

# Clip path and suffixes from filename.
prefix=$(basename $fqgz)
prefix=${prefix/.sorted/}
prefix=${prefix/.fq.gz/}     # That's what the previous step of pipeline uses as suffix
prefix=${prefix/.fastq.gz/}  # Cover this base too.
prefix=${prefix/.trimmed/}
#############



## Module Loading ##
# module load bowtie/1.2.2-foss-2018b
# module load samtools/1.9-foss-2018b
# module load umi-tools/1.0.0-foss-2018b-python-3.6.6
#####################

# Create destination.
if [ ! -d "$outDir" ]; then
	echo "${prefix}: Creating ${outDir}"
  mkdir -p ${outDir}
fi

# Temporary dir
mkdir -p ${outDir}/tmp/${prefix}

echo ""
if [ "$ntnum" -gt 0 ]; then
    if [ "$clipped" -eq 1 ]; then
        # Put UMI back to the front of the read.
        echo "${prefix}: Unclipping UMIs"
        # Grab UMI from end of header.
        # Prepend UMI to sequence.
        # Pad quality string with correct amount of dummy good qualities.
        zcat $fqgz |  perl -e '$hold="";
                              $part=0;
                              while($line = <STDIN>){
                                $part++;
                                if ($part == 1){
                                  $hold=$1 if ( $line =~/([ATGC]+)$/ );
                                  print STDOUT $line;
                                }
                                elsif ($part == 2){
                                  print STDOUT join "", ($hold, $line);
                                }
                                elsif ($part == 3){
                                  print STDOUT $line;
                                }
                                elsif ($part == 4){
                                  print STDOUT join "", ("F" x length($hold), $line);
                                  $part=0;
                                  $hold="";
                                }
                              }' > ${outDir}/tmp/${prefix}/${prefix}_unclipped.fastq
        gzip -f ${outDir}/tmp/${prefix}/${prefix}_unclipped.fastq
        fqgz=${outDir}/tmp/${prefix}/${prefix}_unclipped.fastq.gz
    fi

    echo ""
    echo "${prefix}: Retrieving UMIs"
    repl() {        # Create a UMI pattern of appropriate length.
        printf "N"'%.s' $(seq 1 $1);
    }
    umi_tools extract --temp-dir=${outDir}/tmp/${prefix} -I $fqgz -S ${outDir}/${prefix}_umi-clipped.fastq.gz -p $(repl $ntnum) --extract-method=string --quality-encoding=phred33
    fqgz="$outDir/${prefix}_umi-clipped.fastq.gz"

    if [ "$clipped" -eq 1 ]; then
        rm ${outDir}/tmp/${prefix}/${prefix}_unclipped.fastq.gz
    fi
else
    echo "${prefix}: Skipping UMI extraction"
fi

echo ""
echo "${prefix}: Decompressing FastQ"
gunzip -f -c $fqgz > ${outDir}/${prefix}.fastq

cnt=$(echo "$(wc -l <${outDir}/${prefix}.fastq)/4" | bc)
echo "${prefix} STATS: $cnt unaligned reads"
printf '%s\t%d\n' "unaligned" $cnt >> $metrics


echo ""
echo "${prefix}: Aligning to ${INDEX}"
bowtie --strata -a -m $maxmultis --phred33-quals --chunkmbs 512 -v ${mm_tol} -p 4 -q -S $INDEX ${outDir}/${prefix}.fastq ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sam
samtools view -@ 3 -b ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sam > ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.bam
aln="${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.bam"
rm ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sam

cnt=$(samtools view -c -f 4 $aln)
echo "${prefix} STATS: $cnt reads failed alignment"
printf '%s\t%d\n' "_failed" $cnt >> $metrics

# cnt=$(samtools view -c -F 4 $aln)
${XDIR}/CountMappedReadsBAM.py -b $aln > ${outDir}/${prefix}.cnt
cnt=$(cut -f 3 ${outDir}/${prefix}.cnt)   # Bowtie does not mark secondary alignments, so samtools count is inflated by multimappers
cntm=$(cut -f 5 ${outDir}/${prefix}.cnt)
echo "${prefix} STATS: $cnt aligned reads of which ${cntm} multimappers"
printf '%s\t%d\n' "_aligned" $cnt >> $metrics
printf '%s\t%d\n' "__alignedmutli" $cntm >> $metrics
rm ${outDir}/${prefix}.cnt


echo ""
if [ "$ntnum" -gt 0 ] && [ "$nodedup" -eq 0 ]; then
    echo "${prefix}: Sorting and indexing the aligned BAM"
    samtools sort -@ 3 -m 9800M -o ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.sorted.bam $aln
    rm $aln
    aln="${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.sorted.bam"
    samtools index -@ 3 $aln

    echo ""
    echo "${prefix}: Removing duplicates"
    umi_tools dedup --temp-dir=${outDir}/tmp/${prefix} -I $aln -S ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam --edit-distance-threshold=${mm_umi} --no-sort-output
    aln="${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam"
else
  echo "${prefix}: Skipping duplicate removal"
fi

# cnt=$(samtools view -c -F 4 $aln)
cnt=$(${XDIR}/CountMappedReadsBAM.py -b $aln | cut -f 3)
echo "${prefix} STATS: $cnt deduplicated aligned reads"
printf '%s\t%d\n' "deduplicated" $cnt >> $metrics


rm ${outDir}/${prefix}.fastq
rm -r ${outDir}/tmp/${prefix}  # Leave tmp/, in case there are other instances that also need access to it (processing other samples in parallel).

echo ""
echo "${prefix}: Finished alignment and PCR de-duplication"
exit $?
#############
