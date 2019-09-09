#!/usr/bin/env sh

set -e

source ~/miniconda3/bin/activate mybasics
#module load samtools/1.9-foss-2017a
#module load pysam/0.15.1-foss-2017a-python-3.6.4

in=$1
out=$2
pattern='TTCCAGCATAGCTCTTAAAC'

sequtilities.py T $in --samPatternStats $pattern ${#pattern} N \-4 4 1000000 > $out
bam_barcode_report.R $out $(dirname $out)
