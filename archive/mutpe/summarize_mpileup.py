#!/usr/bin/env python3

from __future__ import print_function

import sys, re, argparse
from collections import Counter

 # Info
usage = "MutPE mutation quantification"
version = "1.3"


parser = argparse.ArgumentParser()
parser.add_argument('-p','--pilein', type=str, help=" Pileup input.")
parser.add_argument('-m','--mutout', type=str, help=" Point mutation stats output.")
parser.add_argument('-i','--indelout', type=str, help=" Insertion & deletion stats output.")
params = parser.parse_args()

pileupFile = params.pilein

mutDict = dict()    # substitutions tallies and deletions
inDict = dict()  # insertion and deletion tallies
coverage = dict()
noise = re.compile(r"\^\S|\$")      # Read starts with mapping quality character, or read ends.
indel = re.compile(r"([-+])(\d+)")   # indel length declarations, capture the type and size of the event.


with open(params.pilein) as f:
    for line in f:
        # Break up the line
        fields = line.rstrip().split("\t")

        cov = fields[3]
        if int(cov) == 0:   # ignore uncovered positions
            continue
        chr = fields[0]
        pos = fields[1]     # samtools mpileup coordinate is 1-based
        ref = fields[2]
        pile = fields[4]
        idpos = str(int(pos) + 0.5)
        if not chr in mutDict:
            mutDict[chr] = dict()
        if not chr in inDict:
            inDict[chr] = dict()
        if not pos in mutDict[chr]:
            mutDict[chr][pos] = Counter()
        if not pos in inDict[chr]:
            inDict[chr][pos] = Counter()
            inDict[chr][idpos] = Counter()
        coverage[pos] = cov
        coverage[idpos] = cov

        # NOTES about the pileup format
        # Insertions : The entire insertion is shown at the position as the preceding match/mismatch.
        #              At these positions a read will have two events happening simulataneously. A naive tally of mutation rates can result in rates up to 200%, causing confusion.
        #              To circumvent that, I print out the insertion stats as decimal positions between the annotated position and the next one. This is OK for my downstream viz
        #              code treating position as a continuous or categorical variable, but will certainly break genome browsers etc tools that expect integer postions.
        # Deletions : Remember that read positions go vertically, not horizontally. A deletion length will be declared at one line, followed by the missing bases,
        #             but the actual deleted bases will be in the following LINES as * characters.

        # Delete read-end symbols, mapping qualities and deletion extensions.
        pile = noise.sub('', pile)
        # Don't care about strands. (presuming the substitutions is given relative to the reference base, not relative to the matching strand)
        pile = pile.upper().replace(',', '.')
        ref = ref.upper()

        # Count deleted bases.
        inDict[chr][pos].update({'-' : pile.count('*')})
        # Remove deletions
        pile = pile.replace('*', '')

        # Count indel declarations by length. Then remove them.
        m = indel.search(pile)
        while(m):
            # I offset declarations by half a position, as:
            #   insertions happen between reference positions anyway and
            #   deletion declarations in the pileup are mapped on the last non-deleted base, which can cause confusion/mistakes and clashes with the actual deleted bases.
            inDict[chr][idpos].update( [m.group(1) + m.group(2)] )
            # Remove insertion and deletion declaration
            pile = re.sub('\\' + m.group(1) + m.group(2) + '[AGCTN]' * int(m.group(2)), '', pile, 1)
            m = indel.search(pile)

        # Sanity check after editing the bases string.
        if not len(pile) <= int(cov):
            sys.stderr.write("Something unexpected may be going on in " + params.pilein + " at pos:\t" + pos + "\tlen:" + str(len(pile)) + "\tcov:" + str(cov) + "\n")

        # Tally up the matches/mismatches [ACGTN.] Mismatches in the ref>mut format. For matches, just the ref base.
        mutDict[chr][pos].update( [ref + '>' + x if x != '.' else ref for x in list(pile)] )

with open(params.mutout, 'w') as mo:
    mo.write("chr\tpos\ttype\tcount\tdepth\n")
    amplicons = sorted(mutDict.keys(), key=str)
    for chr in amplicons:
        positions = sorted(mutDict[chr].keys(), key=int)
        for pos in positions:
            for mut, cnt in mutDict[chr][pos].most_common():
                mo.write(chr + "\t" + pos + "\t" + mut + "\t" + str(cnt) + "\t" + str(coverage[pos]) + "\n")

with open(params.indelout, 'w') as ido:
    ido.write("chr\tpos\ttype\tcount\tdepth\n")
    amplicons = sorted(inDict.keys(), key=str)
    for chr in amplicons:
        positions = sorted(inDict[chr].keys(), key=float)
        for pos in positions:
            for mut, cnt in inDict[chr][pos].most_common():
                ido.write(chr + "\t" + pos + "\t" + mut + "\t" + str(cnt) + "\t" + str(coverage[pos]) + "\n")
