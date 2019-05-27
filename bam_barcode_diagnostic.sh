#!/usr/bin/env sh

source ~/miniconda3/bin/activate mybasics
module load samtools/1.9-foss-2017a

in=$1
out=$2

samtools view $in | sequtilities.py D --samPatternStats TTCCAGCATAGCTCTTAAAC \-4 4 > ${out}

