#!/usr/bin/env python3

import sys
from Bio import SeqIO
from argparse import ArgumentParser, RawDescriptionHelpFormatter

usage = "Chromosome surgery: Splice something into and/or out of a chromosome."
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-f", "--fasta", type=str, required=True, help="Input FASTA.")
parser.add_argument("-o", "--output", type=str, required=True, help="Output FASTA.")
parser.add_argument("-c", "--cid", type=str, required=True, help="Chromosome ID to edit.")
parser.add_argument("-i", "--insert", type=str, required=False, help="FASTA of sequence to insert.")
parser.add_argument("-I", "--incision", type=int, required=False, help="1-based nucleotide after which to insert the insert.")
parser.add_argument("-e", "--excision_start", type=int, required=False, help="1-based nucleotide that is the first to delete (0).")
parser.add_argument("-E", "--excision_end", type=int, required=False, help="1-based nucleotide that is the last to delete (0).")
args = parser.parse_args()


# Harmless defaults
splice_in = ''


# Get insert
if args.insert:
    with open(args.insert, 'r') as splicein:
        record = list(SeqIO.parse(splicein, 'fasta'))[0]
        splice_in = record.seq
    
    # No need to shift the incision coordinate. 
    # The 1-based right-closed index after which to cut is the same location as the 0-based right-open substring end before the cut.

if args.excision_start and args.excision_end:
    excision_start = args.excision_start
    excision_end = args.excision_end
    # Pythonize start coordinate from 1-based left-closed to 0-based left-closed.
    excision_start -= 1
    # No need to change the end coordinate. The 1-based right-closed index is the same location as the 0-based right-open substring end.

if args.insert and args.excision_start:
    # Do excision after the incision. 
    # Adjust coordinates.
    if args.excision_start > args.incision and args.excision_end > args.incision:
        excision_start = args.excision_start + len(splice_in)
        excision_end = args.excision_end + len(splice_in)
   
    elif args.excision_start < args.incision and args.excision_end < args.incision:
        pass    # The incision will be applied first, no need to adjust it. The excision is unaffected by the incision anyway.
    else:
        sys.err.write('Error: Cannot apply the specified coordinates. Excision end must be after excision start, and the incision cannot be inside the excision.')
    

# Parse and apply edit
with open(args.fasta, 'r') as genome,   open(args.output, 'w') as out:
    for record in SeqIO.parse(genome, 'fasta'):

        # Only edit the relevant entry
        if (record.id == args.cid):
            # Splice-in
            record.seq = record.seq[:args.incision] + splice_in + record.seq[args.incision:]
            # Splice-out
            record.seq = record.seq[:excision_start] + record.seq[excision_end:]

        # Output all the entries    
        SeqIO.write(record, out, 'fasta')

print("Done")
