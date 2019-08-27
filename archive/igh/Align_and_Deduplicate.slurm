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
# It takes as input a SINGLE-ENDED fastq file, the name of the dir where output files will be created,
# and a precomputed bowtie index.
#
# Last reviewed: 05/aug/2019	by: kimon.froussios@imp.ac.at
####################

set -e

## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -f FASTQ_FILE -o OUTDIR -i INDEX -r CHR_SIZES -l LOCUS [-m MISMATCHES] [-n UMI_LENGTH] [-u CLIPPED]"
    exit 1
}
# Defaults.
mm_tol=1
ntnum=10
mm_umi=1
clipped=0
# Parse options.
while getopts 'f:o:i:l:r:m:n:u:' flag; do
  case "${flag}" in
    f) fqgz="${OPTARG}" ;;			# Compressed FASTQ file.
    o) outDir="${OPTARG}" ;;    # Output directory.
    i) INDEX="${OPTARG}" ;;			# Bowtie index prefix. (not bowtie2)
    r) ref_sizes="${OPTARGS}" ;;  # Table of reference sequence sizes
    l) LOCUS="${OPTARG}" ;;			# Short nickname for the index used, to be integrated into output file names.
    m) mm_tol="${OPTARG}" ;;		# Mismatch tolerance. Choose between 0 and 3 (1).
    n) ntnum="${OPTARG}" ;;			# Length of UMI at the start of reads (10). Use 0 to skip deduplication.
    u) clipped="${OPTARG}" ;;     # 1: UMI already clipped and appended to end of read title, 0: UMI not clipped
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
#############



## Module Loading ##
module load bowtie/1.2.2_p1-foss-2017a
module load samtools/1.9-foss-2017a
module load umi-tools/1.0.0-foss-2017a-python-3.6.4
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
  # Put UMI back to the front of the read.
  if [ "$clipped" -eq 1 ]; then
    echo "${prefix}: Unclipping UMIs"
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
else
  echo "${prefix}: Skipping UMI extraction"
fi

echo ""
echo "${prefix}: Decompressing FastQ"
gunzip -f -c $fqgz > ${outDir}/${prefix}.fastq

echo "${prefix}: Aligning to ${INDEX}"
bowtie --best --strata -m 60 --phred33-quals --chunkmbs 256 -v ${mm_tol} -a -p 3 -q -S $INDEX ${outDir}/${prefix}.fastq ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sam

echo "${prefix}: Sorting and indexing the aligned BAM"
samtools view -@ 3 -b ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sam > ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.bam
samtools sort -@ 3 -m 9800M -o ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sorted.bam ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.bam
samtools index -@ 3 ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sorted.bam

echo ""
if [ "$ntnum" -gt 0 ]; then
  echo "${prefix}: Removing duplicates"
  umi_tools dedup --temp-dir=${outDir}/tmp/${prefix} -I ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sorted.bam -S ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam --edit-distance-threshold=${mm_umi} --no-sort-output
  
  echo "${prefix}: Sorting and indexing the deduplicated BAM"
  samtools sort -@ 3 -m 9800M -o ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.deduped.sorted.bam ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam
  samtools index -@ 3 ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.deduped.sorted.bam
else
  echo "${prefix}: Skipping duplicate removal"
fi

rm ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.sam
rm ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.aln.bam*
if [ "$ntnum" -gt 0 ]; then
  rm ${outDir}/${prefix}_${LOCUS}_mm${mm_tol}.deduped.bam*
fi
rm ${outDir}/${prefix}.fastq
rm -r ${outDir}/tmp/${prefix}  # Leave tmp/, in case there are other instances that also need access to it (processing other samples in parallel).
echo "${prefix}: Finished alignment and PCR de-duplication"

module unload bowtie
module unload samtools
module unload umi-tools

exit $?
#############
