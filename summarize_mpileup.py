#!/usr/bin/env python3

from __future__ import print_function

import sys, re
from collections import Counter

 # Info
usage = "MutPE mutation quantification"
version = "1.2"

pileupFile = sys.argv[1]

mutDict = dict()    # substitutions tallies and deletions
inDict = dict()  # insertion tallies
coverage = dict()
noise = re.compile(r"\^\S|\$|\*")      # Read starts with mapping quality character, read ends, and deletion extensions.
insertions = re.compile(r"[+](\d+)")   # insertion length declarations
deletions = re.compile(r"-(\d+)")      # deletion length declarations

with open(pileupFile) as f:
    for line in f:
        # Break up the line
        fields = line.rstrip().split("\t")
        
        cov = fields[3]
        if int(cov) == 0:   # ignore uncovered positions
            continue
        chr = fields[0]
        pos = fields[1]
        ref = fields[2]
        pile = fields[4]
        if not chr in mutDict:
            mutDict[chr] = dict()
        if not chr in inDict:
            inDict[chr] = dict()
        if not pos in mutDict[chr]:
            mutDict[chr][pos] = Counter()  # It might have been created in the previous line, because of how deletions are defined in pileup.
        coverage[pos] = cov
        
        # NOTES about the pileup format
        # Insertions : The entire insertion is shown at the position as the preceding match/mismatch. 
        #              At these positions a read will have two events happening simulataneously. A naive tally of mutation rates can result in rates up to 200%, causing confusion.
        #              To circumvent that, I print out the insertion stats as decimal positions between the annotated position and the next one. This is OK for my downstream viz 
        #              code treating position as a continuous numeric or a categorical variable, but will certainly break genome browsers etc tools that expect integer postions.
        # Deletions : Remember that read positions go vertically, not horizontally. A deletion length will be declared at one line, followed by the missing bases, 
        #             but the actual deleted bases will be in the following LINES as * characters.
        #             I could assign the deletion to the first deleted base, but I think I'l
        
        # Delete read-end symbols, mapping qualities and deletion extensions.
        pile = noise.sub('', pile)
        
        # Find insertions and reduce them to the leading + character, discarding length and inserted sequence.
        for m in insertions.finditer(pile):
            pile = re.sub(m.group(1) + '.{' + m.group(1) + '}', '', pile, 1)
        # Isolate the insertions
        pile2 = re.sub(r"[^+]", '', pile)
        # Remove insertions
        pile = re.sub(r"[+]", '', pile)
        
        # Find deletion declarations and reduce them to the leading - character, discarding length and sequence.
        for m in deletions.finditer(pile):
            pile = re.sub(m.group(1) + '.{' + m.group(1) + '}', '', pile, 1)
        # Isolate the deletion declarations.
        pile3 = re.sub(r"[^\-]", '', pile)
        # Remove deletion declarations.
        pile = re.sub(r"[\-]", '', pile)
        
        # Sanity check after editing the bases string.
        if not len(pile) <= int(cov):
            sys.stderr.write("Something unexpected may be going on at pos:\t" + pos + "\tlen:" + str(len(pile)) + "\tcov:" + str(cov) + "\n")
        
        # Don't care about strands
        pile = pile.upper().replace(',', '.')
        # Tally up the matches/mismatches [ACGTN.]
        mutDict[chr][pos].update(pile)
        # Tally up the insertions [+]
        if len(pile2) > 0:
            inDict[chr][pos] = Counter()
            inDict[chr][pos].update(pile2)
        # Tally up the deletions [-]
        if len(pile3) > 0:
            if not int(pos) + 1 in mutDict[chr]:
                mutDict[chr][str(int(pos) + 1)] = Counter()
                mutDict[chr][str(int(pos) + 1)].update(pile3)

# Output header
print("seq\tpos\ttype\tcount\tdepth")            

amplicons = sorted(mutDict.keys(), key=str)

for chr in amplicons:
    positions = sorted(mutDict[chr].keys(), key=int)
    for pos in positions:
        for mut, cnt in mutDict[chr][pos].most_common():
            print(chr + "\t" + pos + "\t" + mut + "\t" + str(cnt) + "\t" + str(coverage[pos]))
        if pos in inDict[chr]:
            for mut, cnt in inDict[chr][pos].most_common():      # in case I change my mind and log other things there too, or if I log insertion lengths or sequence, instead of just the presence.
                # Print at a decimal offset position, to prevent overlapping the base that occupies the actual position.
                print(chr + "\t" + str(int(pos) + 0.5) + "\t" + mut + "\t" + str(cnt) + "\t" + str(coverage[pos]))
