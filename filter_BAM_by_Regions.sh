#!/usr/bin/sh

bam=$1      # input BAM
regions=$2  # SAM header file, with only the sequences to keep
out=$3      # output BAM

# module load samtools/1.9-foss-2017a

samtools view $bam | sequtilities.py D --samFltrRegs $regions | samtools view -b - > $out
