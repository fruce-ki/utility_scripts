from Bio import SeqIO



from argparse import ArgumentParser, RawDescriptionHelpFormatter

usage = "Chromosome surgery: Splice something into and/or out of a chromosome."
# Main Parsers
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-f", "--fasta", type=str, required=True, help="Genomic FASTA.")
parser.add_argument("-o", "--cid", type=str, required=True, help="Chromosome ID to edit.")
parser.add_argument("-i", "--insert", type=str, required=False, help="FASTA of sequence to insert.")
parser.add_argument("-I", "--incision", type=str, required=False, help="1-based nucleotide after which to insert the insert.")
parser.add_argument("-e", "--excision_start", type=str, required=False, help="1-based nucleotide that is the first to delete.")
parser.add_argument("-E", "--excsiion_end", type=str, required=False, help="1-based nucleotide that is the last to delete.")
args = parser.parse_args()

splicein = ''

if (args.insert):
    with open(args.insert, 'r') as injection:
        SeqIO.parse(injection, 'fasta')


with open(args.fasta, 'r') as genome:
    for record in SeqIO.parse(genome, 'fasta'):
        if (record.id == args.cid):

            pass