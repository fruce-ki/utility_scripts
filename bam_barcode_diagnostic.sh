#!/usr/bin/env sh

set -e

source ~/miniconda3/bin/activate mybasics
module load samtools/1.9-foss-2017a

in=$1
out=$2
pattern='TTCCAGCATAGCTCTTAAAC'

samtools view $in | sequtilities.py D --samPatternStats $pattern ${#pattern} N \-4 4 > ${out}
bam_barcode_report.R $out
