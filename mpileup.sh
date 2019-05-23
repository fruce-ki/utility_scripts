#!/usr/bin/env sh

module load samtools/1.9-foss-2017a

fa=$1 ; shift
bam=$1 ; shift
params=( "$@" )

samtools mpileup -A -B -f $fa -d 1000000000 ${bam} > ${bam/.bam/.pileup}
