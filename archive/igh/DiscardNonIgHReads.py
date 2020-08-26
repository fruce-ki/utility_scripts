#!/usr/bin/env python3

import re, os, sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import Counter

import pysam
from intervaltree import Interval, IntervalTree
from pybedtools import BedTool

#############
## Read Me ##
#############
# Complete redesign of the discard step by Tobias Neumann.
#
# 11/10/2019 Extensively reworked for simpler logic by Kimon Froussios.
#
# The major assumption of this script is that alignment was done with bowtie --best --strata -a
# This means that only equally-best alignments are reported as multimappers (by number of mismatches).
#
# Last revised: 26/jun/2020   by: kimon.froussios@imp.ac.at

#################
## Parameters ##
################

usage = "Filter script to retain IgH multimappers"
version = "0.3.0"
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-a", "--allSam", type=str, required=True, dest="allSamFile", help="BAM file with all the mapped reads")
parser.add_argument("-g", "--genome", type=str, required=True, dest="genome", help="The chromosome name. Also serves to automatically configure relevant settings.")
parser.add_argument("-o", "--outDir", type=str, required=True, dest="outDir", help="Destination directory.")
args = parser.parse_args()

#################
## Pre-defined Settings
################

IgHTree = dict()
annotfile = None
if args.genome == "B18_DFR" or args.genome == "CH12_VDJ":
    IgHChr = "12"
    IgHTree[114496979:117248166] = 0
    #annotfile="/groups/pavri/bioinfo/Marisol/custom_seqs/B18_DFR_newest/filter/IgH_repeat_regions_withPromoter_noChr.bed"
    annotfile="/groups/pavri/bioinfo/Marisol/custom_seqs/B18_DFR_newest/filter/refseq.mm9.igh.merged.noChr.bed"

elif args.genome == "Ramos_IgH":
    IgHChr = ["chr14"]
    #IgHTree[106032614:107288051] = 0
    IgHTree["chr14"] =IntervalTree()
    IgHTree["chr14"][105836764:106875071] = 0
    #annotfile="/groups/pavri/bioinfo/Marisol/custom_seqs/Ramos_IgH_newest/filter/Ramos_IGH_filter_regions.bed"
    # annotfile="/groups/zuber/zubarchive/USERS/tobias/rushad/custom_seqs/Ramos_IgH_hg38/filter/hg38_vdj_exons.bed"
    annotfile="/users/kimon.froussios/pavri/ursi/Ramos_PRO/aux/Hg38_plus_RamosIgH/hg38_vdj_exons.bed"

elif args.genome in ['B1-8hi_mm10_190516_181014-most-recent',
                     'B18hi_HDRc_BfaI_ctrl_190513_US',
                     'B18hi_HDR1_dTATA_190513_US',
                     'B18hi_HDR2_dTATA_+1T_190513_US',
                     'B18hi_HDR3_TATTAC_190513_US',
                     'B18hi_HDR4_TATTAC_+1T_190513_US',
                     'CH12_VDJ_200122']:
    IgHChr = ["NC_000078.6", "NT_114985.3"]   # NCBI mm10 chr12 and unplaced alternative scaffold
    IgHTree["NC_000078.6"] = IntervalTree()
    IgHTree["NC_000078.6"][113572929:116009954] = 0  # mm10 main V regions (183)
    IgHTree["NC_000078.6"][113407535:113526809] = 0  # mm10 D regions (21)
    IgHTree["NC_000078.6"][113428514:113429833] = 0  # mm10 J regions (4)
    IgHTree["NT_114985.3"] = IntervalTree()
    IgHTree["NT_114985.3"][466661:1711874] = 0      # mm10 alternate scaffold V region (11)
    IgHTree["NT_114985.3"][277878:300620] = 0      # mm10 alternate scaffold D region (1)
    annotfile = "/groups/pavri/Kimon/ursi/PRO/aux/mm10_ncbi/GCF_000001635.26_GRCm38.p6_VDJ.bed"

else:
    raise Exception("The locus reference name " + args.genome + " is not among the pre-defined ones!")


##########
# Categorize read names
##########

LOCUS = Counter()  # might map multiple times on the artificial chromosome.
IGH = Counter()    # the locus was built form IgH sequences. So if reads map exactly once to IgH, I can include themas locus. If more than once, then too generic.
other = set()      # mere existance of a match outside IgH is enough to disqualify a read.

inSam = pysam.AlignmentFile(args.allSamFile, "rb")
# Multiple alignments of a read are reported individually.
for read in inSam:
    if not read.is_unmapped:
        if read.reference_name == args.genome:
            # Mapped to target locus.
            LOCUS.update([read.query_name])
        elif read.reference_name in IgHChr and len(IgHTree[read.reference_name][read.reference_start:read.reference_end]) > 0:
            # Overlaps native IgH regions
            IGH.update([read.query_name])
        else:
            # Anywhere else
            other.add(read.query_name)
inSam.close()

# Filter reads
# unique : maps to locus once and up to once in IgH (the V,D,J from which the artificial chromosome was assembled), but nowhere outside IgH.
# multi : maps to locus once and more than once in IgH (ie. more generic across Vs Ds and Js), but nowhere outside IgH.
unique = [x for x in LOCUS if LOCUS[x] == 1 and (x not in IGH or IGH[x] == 1) and x not in other]
multi = [x for x in LOCUS if LOCUS[x] == 1 and x in IGH and IGH[x] > 1 and x not in other]

##########
# Categorize reads
##########

prefix = re.sub(".sam|.bam", "", os.path.basename(args.allSamFile))
prefix = re.sub(".deduped|.aln|_subject", "", prefix)

# Re-open input from the start.
inSam = pysam.AlignmentFile(args.allSamFile, "rb")
# Mapped only once to target locus and up to once in genomic IgH.
outUnique = pysam.AlignmentFile(os.path.join(args.outDir, prefix + "_unique.bam"), "wb", template=inSam)
# Mapped once or more to target locus and more than once to genomic IgH.
outNonUnique = pysam.AlignmentFile(os.path.join(args.outDir, prefix + "_nonunique.bam"), "wb", template=inSam)
# Unique and non-unique mapped to locus (ie. the above two together).
outLocus = pysam.AlignmentFile(os.path.join(args.outDir, prefix + "_locus.bam"), "wb", template=inSam)

for read in inSam:
    # Discard unmapped reads.
    if not read.is_unmapped:
        if (read.reference_name == args.genome):
            if read.query_name in unique:
                outUnique.write(read)
                outLocus.write(read)
            elif read.query_name in multi:
                outNonUnique.write(read)
                outLocus.write(read)

inSam.close()
outUnique.close()
outNonUnique.close()
outLocus.close()

