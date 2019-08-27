#!/usr/bin/env python2.7

from __future__ import print_function
import sys, os

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import commands, re, os, sys, time,subprocess,textwrap
import pysam, random
from collections import defaultdict

from intervaltree import Interval, IntervalTree
from pybedtools import BedTool

#############
## Read Me ##
#############
# Complete redesign of the discard step by Tobias Neumann.

# Updated to fit the redesigned workflow by Kimon Froussios.

#################
## Parameters ##
################

def bedToIntervallTree(bed):

    tree = {}

    exons = BedTool(bed)

    for exon in exons:

        if (not tree.has_key(exon.chrom)) :
            tree[exon.chrom] = IntervalTree()

        tree[exon.chrom][exon.start:(exon.end + 1)] = 0

    return tree

usage = "Filter script to retain IgH multimappers"
version = "0.2.0"
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter, version=version)

parser.add_argument("-c", "--chromosomeSam", type=str, required=True, dest="chromosomeFile", help="SAM file with reads from the custom chromosome only")
parser.add_argument("-a", "--allSam", type=str, required=True, dest="allSamFile", help="BAM file with all the mapped reads")
parser.add_argument("-g", "--genome", type=str, required=True, dest="genome", help="The name of the custom chromosome.")

args = parser.parse_args()

IgHTree = IntervalTree()

annotfile = None
if args.genome == "B18_DFR" or args.genome == "CH12_VDJ":
    IgHChr = "12"
    IgHTree[114496979:117248166] = 0
    #annotfile="/groups/pavri/bioinfo/Marisol/custom_seqs/B18_DFR_newest/filter/IgH_repeat_regions_withPromoter_noChr.bed"
    annotfile="/groups/pavri/bioinfo/Marisol/custom_seqs/B18_DFR_newest/filter/refseq.mm9.igh.merged.noChr.bed"

elif args.genome == "Ramos_IgH":
    IgHChr = "chr14"
    #IgHTree[106032614:107288051] = 0
    IgHTree[105836764:106875071] = 0
    #annotfile="/groups/pavri/bioinfo/Marisol/custom_seqs/Ramos_IgH_newest/filter/Ramos_IGH_filter_regions.bed"
    annotfile="/groups/zuber/zubarchive/USERS/tobias/rushad/custom_seqs/Ramos_IgH_hg38/filter/hg38_vdj_exons.bed"

elif args.genome in ['B1-8hi_mm10_190516_181014-most-recent', 
                     'B18hi_HDRc_BfaI_ctrl_190513_US', 
                     'B18hi_HDR1_dTATA_190513_US', 
                     'B18hi_HDR2_dTATA_+1T_190513_US', 
                     'B18hi_HDR3_TATTAC_190513_US', 
                     'B18hi_HDR4_TATTAC_+1T_190513_US']:
    IgHChr = "NC_000078.6"            # NCBI mm10 chr12
    IgHTree[113572929:116009954] = 0  # mm10 main V regions
    #IgHTree[466661:1711874] = 0      # mm10 alternate scaffold V region
    IgHTree[113430528:113526809] = 0  # mm10 D regions
    IgHTree[113428514:113429833] = 0  # mm10 J regions
    annotfile = "/groups/pavri/Kimon/ursi/PRO/aux/mm10_ncbi/GCF_000001635.26_GRCm38.p6_VDJ.bed"

else:
    raise Exception("The locus reference name is not among the pre-defined ones!")

# Dics to keep reads

chrDict = {}
singletonDict = {}
uniqueDict = {}
exonDict = {}
IgHDict = {}
discardedDict = {}

# Read filter regions to IntervalTree

IgHExons = bedToIntervallTree(annotfile)

# Read reads mapping to extra-chromosome

chrSam = pysam.AlignmentFile(args.chromosomeFile, "rb")
for read in chrSam:
    if (not read.is_unmapped and chrSam.getrname(read.reference_id) == args.genome) :
        # if it is already already seen, remove it from the singletons
        if (singletonDict.has_key(read.query_name)) :
            singletonDict.pop(read.query_name,None)
        # not seen before, add it as singleton
        else :
            singletonDict[read.query_name] = read
        # in any case, mark it as mapping to the artificial locus
        chrDict[read.query_name] = read

# Open all different diagnostic sams

outfileUnique = pysam.AlignmentFile(re.sub(".sam","",args.chromosomeFile)+"_unique.bam", "wb", template=chrSam)
outfileExon = pysam.AlignmentFile(re.sub(".sam","",args.chromosomeFile)+"_exonfiltered.bam", "wb", template=chrSam)
outfileIgH = pysam.AlignmentFile(re.sub(".sam","",args.chromosomeFile)+"_ighfiltered.bam", "wb", template=chrSam)
outfileDiscarded = pysam.AlignmentFile(re.sub(".sam","",args.chromosomeFile)+"_discarded.bam", "wb", template=chrSam)
outfileMultiple = pysam.AlignmentFile(re.sub(".sam","",args.chromosomeFile)+"_multiple.bam", "wb", template=chrSam)
outfileNorealignment = pysam.AlignmentFile(re.sub(".sam","",args.chromosomeFile)+"_noRealign.bam", "wb", template=chrSam)
outfileAdded = pysam.AlignmentFile(re.sub(".sam","",args.chromosomeFile)+"_added.bam", "wb", template=chrSam)

chrSam.close()

# For before and after shots

readsWithoutRealignment = dict()
readsAdded = dict()

# Walk through genome-wide alignments

allSam = pysam.AlignmentFile(args.allSamFile, "rb")

for read in allSam:
    # Every mapped read that is on the correct chromosome and is not already mapped to the artificial locus
    if (not read.is_unmapped and allSam.getrname(read.reference_id) != args.genome and chrDict.has_key(read.query_name)) :

        chr = allSam.getrname(read.reference_id)

        # Check if read is contained in IgH region
        if IgHExons.has_key(chr) :

            query = IgHExons[chr][read.reference_start:read.reference_end]
            #if (len(query) > 0 and not IgHDict.has_key(read.query_name) and not discardedDict.has_key(read.query_name)) :
            if (len(query) > 0 and not discardedDict.has_key(read.query_name)) :
                if (not read.query_name in readsWithoutRealignment and not read.query_name in readsAdded) :
                    readsWithoutRealignment[read.query_name] = read
                else :
                    readsAdded[read.query_name] = read
                    readsWithoutRealignment.pop(read.query_name,None)

                exonDict[read.query_name] = read
                IgHDict.pop(read.query_name,None)
            else :
                query = IgHTree[read.reference_start:read.reference_end]

                #if (len(query) > 0 and not discardedDict.has_key(read.query_name)) :
                if (len(query) > 0 and not discardedDict.has_key(read.query_name)) :
                    if (not exonDict.has_key(read.query_name)) :
                        IgHDict[read.query_name] = read
                        #    exonDict.pop(read.query_name,None)
                    if (not read.query_name in readsWithoutRealignment and not read.query_name in readsAdded) :
                        readsWithoutRealignment[read.query_name] = read
                    else :
                        readsAdded[read.query_name] = read
                        readsWithoutRealignment.pop(read.query_name,None)
                else :
                    discardedDict[read.query_name] = read

                    #chrDict.pop(read.query_name,None)
                    exonDict.pop(read.query_name,None)
                    IgHDict.pop(read.query_name,None)
                    readsAdded.pop(read.query_name,None)
                    readsWithoutRealignment.pop(read.query_name,None)

        else :
            discardedDict[read.query_name] = read
            #chrDict.pop(read.query_name,None)
            exonDict.pop(read.query_name,None)
            IgHDict.pop(read.query_name,None)
            readsAdded.pop(read.query_name,None)
            readsWithoutRealignment.pop(read.query_name,None)

for readname in singletonDict.keys():

    if (readname in readsWithoutRealignment) :
        if (readname in readsAdded) :
            print("OMG read " + readname + " found in multiple dicts!",file=sys.stderr)
            sys.exit(-1)
        outfileNorealignment.write(singletonDict[readname])

    if (readname in readsAdded) :
        if (readname in readsWithoutRealignment) :
            print("OMG read " + readname + " found in multiple dicts!",file=sys.stderr)
            sys.exit(-1)
        outfileAdded.write(singletonDict[readname])

    if (exonDict.has_key(readname)) :
        outfileExon.write(singletonDict[readname])
        if (IgHDict.has_key(readname) or discardedDict.has_key(readname)) :
            print("OMG read " + readname + " found in multiple dicts!",file=sys.stderr)
            sys.exit(-1)
    elif (IgHDict.has_key(readname)) :
        outfileIgH.write(singletonDict[readname])
        if (discardedDict.has_key(readname)) :
            print("OMG read " + readname + " found in discarded dict!",file=sys.stderr)
            sys.exit(-1)

    elif (discardedDict.has_key(readname)) :
        outfileDiscarded.write(singletonDict[readname])

    else :
        outfileUnique.write(singletonDict[readname])
        outfileNorealignment.write(singletonDict[readname])

for readname in chrDict.keys():

    if (exonDict.has_key(readname)) :
        outfileMultiple.write(chrDict[readname])
        if (IgHDict.has_key(readname) or discardedDict.has_key(readname)) :
            print("OMG read " + readname + " found in multiple dicts!",file=sys.stderr)
            sys.exit(-1)
    elif (IgHDict.has_key(readname)) :
        outfileMultiple.write(chrDict[readname])
        if (discardedDict.has_key(readname)) :
            print("OMG read " + readname + " found in discarded dict!",file=sys.stderr)
            sys.exit(-1)

    elif (discardedDict.has_key(readname)) :
        pass

    else :
        outfileMultiple.write(chrDict[readname])

allSam.close()
outfileUnique.close()
outfileExon.close()
outfileIgH.close()
outfileDiscarded.close()
outfileMultiple.close()
outfileNorealignment.close()
outfileAdded.close()

#pysam.sort(re.sub(".sam","",args.chromosomeFile)+"_unique.bam", re.sub(".sam","",args.chromosomeFile)+"_unique_sorted")
#pysam.sort(re.sub(".sam","",args.chromosomeFile)+"_exonfiltered.bam", re.sub(".sam","",args.chromosomeFile)+"_exonfiltered_sorted")
#pysam.sort(re.sub(".sam","",args.chromosomeFile)+"_ighfiltered.bam", re.sub(".sam","",args.chromosomeFile)+"_ighfiltered_sorted")
#pysam.sort(re.sub(".sam","",args.chromosomeFile)+"_discarded.bam", re.sub(".sam","",args.chromosomeFile)+"_discarded_sorted")
