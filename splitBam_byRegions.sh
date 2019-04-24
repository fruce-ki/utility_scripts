#!/usr/bin/sh

bam=$1      # input BAM
regions=$2  # SAM header file, with only the sequences to keep
out=$3      # output BAM

samtools view $bam | sequtilities.py D --samFltrRegs $regions | samtools view -b - > $out
