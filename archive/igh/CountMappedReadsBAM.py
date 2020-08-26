#!/usr/bin/env python3

import re, os, sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import Counter

import pysam

#############
## Read Me ##
#############
# Complete redesign of the discard step by Tobias Neumann.
#
# 22/06/2020 Created
#
# This just scans through the BAM, and reports the number of unique IDs for mapped reads. Unpapped reads will be ignored.
#
# Last revised: 26/jun/2020   by: kimon.froussios@imp.ac.at

#################
## Parameters ##
################

usage = "Count unique read IDs in a BAM"
version = "0.2.0"
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b", "--bam", type=str, required=True, dest="allSamFile", help="BAM file with all the mapped reads")
args = parser.parse_args()

##########
# Collect read names
##########

readID = Counter()

inSam = pysam.AlignmentFile(args.allSamFile, "rb")
# Multiple alignments of a read are reported individually.
for read in inSam:
    if not read.is_unmapped:
        readID.update([read.query_name])
inSam.close()

i = 0
for k in readID.items():
    if k[1] > 1:
        i += 1 

sys.stdout.write(args.allSamFile + '\tReads:\t' + str(len(readID)) + '\tof which multimappers:\t' + str(i) + '\n')
