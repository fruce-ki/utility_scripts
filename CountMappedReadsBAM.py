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

usage = "Count unique read IDs in a BAM. Partially overlapping coordinates are included."
version = "0.3.0"
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b", "--bam", type=str, nargs='+', required=True, dest="allSamFile", help="BAM file with all the mapped reads")
parser.add_argument("-c", "--chr", type=str, required=False, dest="chromosome", help="chromosome")
parser.add_argument("-s", "--start", type=int, required=False, dest="coordstart", help="from coordinate, 0-based, inclusive")
parser.add_argument("-e", "--end", type=int, required=False, dest="coordend", help="to coordinate, 0-based ,not inclusive")
args = parser.parse_args()

##########
# Collect read names
##########

for f in args.allSamFile:
    readID = Counter()

    inSam = pysam.AlignmentFile(f, "rb")
    # Multiple alignments of a read are reported individually.
    for read in inSam:
        if not read.is_unmapped:
            flag = True   # keeping track of if all conditions are individually met, allowing for any combination of conditions to be specified
            if args.chromosome is not None and read.reference_name != args.chromosome:  # reference sequence outside specification
                flag = False
            if args.coordstart is not None and read.reference_end < args.coordstart:    # read mapped before specified location
                flag = False
            if args.coordend is not None and read.reference_start > args.coordend:      # read mapped after specified location
                flag = False        
            
            if flag:
                readID.update([read.query_name])
    inSam.close()

    i = 0
    for k in readID.items():
        if k[1] > 1:
            i += 1 

    sys.stdout.write(f + '\tReads:\t' + str(len(readID)) + '\tof which multimappers:\t' + str(i) + '\n')
