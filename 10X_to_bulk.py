#!/usr/bin/env python3

## Spice 10X R1 (cellcode and UMI) onto R2 (actual 3' read).

import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from argparse import ArgumentParser, RawDescriptionHelpFormatter

usage = "Chromosome surgery: Splice something into and/or out of a chromosome."
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-1", "--R1", type=str, required=True, help="fastq.gz with barcode/umi reads. Will form the 1st part of the merged read.")
parser.add_argument("-2", "--R2", type=str, required=True, help="fastq with sequence reads. Will form the 2nd part of the merged read.")
parser.add_argument("-o", "--out", type=str, required=True, help="Destination fastq.gz for the merged reads.")
args = parser.parse_args()


count=0

with open(args.out, 'w') as out_handle:
    with gzip.open(args.R1, 'rt') as umi_handle:
        with gzip.open(args.R2, 'rt') as seq_handle:

            shi = FastqGeneralIterator(seq_handle)
            ## Parse UMI
            for title1, seq1, qual1 in FastqGeneralIterator(umi_handle):
                ## Parse corresponding sequence. The files are expected to have the same numbeer of entries and in the same order.
                title2, seq2, qual2 = next(shi)

                count += 1

                ## Verify titles match
                if title1.split(' ')[0] != title2.split(' ')[0]:
                    print("Error: Stumbled on unmatched reads at entry %d!" % count)
                    print(title1)
                    print(title2)
                    exit(1)
                else:
                    newtitle = title1.rstrip("\n")
                    newseq = seq1.rstrip("\n") + seq2.rstrip("\n")
                    newqual = qual1.rstrip("\n") + qual2.rstrip("\n")
                    
                    out_handle.write("@%s\n%s\n+\n%s\n" % (newtitle, newseq, newqual))


print("Done!\n%d reads merged." % count)

