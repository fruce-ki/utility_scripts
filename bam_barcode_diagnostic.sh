#!/usr/bin/env sh

source ~/miniconda3/bin/activate mybasics
module load samtools/1.9-foss-2017a

in=$1
out=$2
pattern='TTCCAGCATAGCTCTTAAAC'

samtools view $in | sequtilities.py D --samPatternStats $pattern $(expr ${#pattern} / 2) N \-4 4 0.01 > ${out}
