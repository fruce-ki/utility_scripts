#!/usr/bin/env sh

# split_BAM_by_numMismatches.sh <path/file> <MIN#> <MAX#>

samtools view -h $1 | perl ~/utility_scripts/split_SAM_by_numMismatches.pl $2 $3 | samtools view -b -o ${1/.bam/_$2-$3.bam}
samtools index ${1/.bam/_$2-$3.bam}
