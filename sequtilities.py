#!/usr/bin/env python3

"""sequtilities.py

Author: Kimon Froussios
Date last revised: 04/06/2019

Library of utility functions relevant to sequencing tasks.
"""


import os, sys, argparse, re, csv
import pandas as pd
from collections import Counter

import Levenshtein as lev
import pysam

import fileutilities as fu
import mylogs as ml


def collect_starFinalLogs(flist, all=False):
    """Combine the listed Log.final.out files into a pandas Dataframe.

    File identifiers (filenames) will be trimmed at '_Log'.

    Args:
        flist: A list/FilesList of input files.
        all(bool): Show all fields (False).
    Returns:
        pandas.DataFrame
    """
    rows = None
    if all:
        # Still discarding some irrelevant stuff and blank columns.
        rows = [2,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,27,28,29]
    else:
        rows = [4, 5, 8, 9, 23, 25, 27, 28, 29]
    df = fu.get_crosspoints(flist, cols=[1], rows=rows, colSep=["\|"], header=False, index=0, merge=True)[0]

    spaces = re.compile("\s{2,}|\t")
    quotes = re.compile("\"")
    # Clean up padding from cells.
    for r in range(0, len(df.index)):
        for c in range(0, len(df.columns)):
            df.iloc[r,c] = spaces.sub("", str(df.iloc[r,c]))
    # Clean up field descriptions.
    index = df.index.values.tolist()
    for i in range(0,len(index)):
        index[i] = quotes.sub("", index[i])
        index[i] = spaces.sub("", index[i])
    df.index = index
    # Clean up file identifiers from suffixes.
    columns = df.columns.values.tolist()
    for c in range(0, len(columns)):
        columns[c] = str(columns[c]).split("_Log")[0]
    df.columns = columns
    # Transpose and add name to new index column.
    df = pd.DataFrame.transpose(df)
    df.index.name = "Name"
    return df


# Helper function
def gtf2pandas(flist):
    """Import GTF files as pandas.DataFrames.

    Two columns will be added: the parent_id (gene) and target_id (transcript),
    as extracted from the attributes column. The attributes column itself will remain as is.
    The choice of column names was made for compatibility with sleuth.

    Args:
        flist: A list/FilesList of GTF files.
    Returns:
        [pandas.DataFrame]: List of dataframes, one per file.
    """
    # Use my own parser, because it already handles input from multiple files or STDIN.
    input = fu.get_columns(flist, cols=list(range(0,9)), colSep=["\t"], header=False, index=None, merge=False)
    result = []
    for gtf in input:
        gtf.columns = ["chr", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"]
        gtf["parent_id"] = gtf["attributes"].str.extract('gene_id \"?([^\";]+)', expand=False)
        gtf["target_id"] = gtf["attributes"].str.extract('transcript_id \"?([^\";]+)', expand=False)
        gtf.set_index(gtf["parent_id"], inplace=True)
        result.append(gtf)
    return result

def gtf2premrna(gtfs, filter=True):
    """Infer pre-mRNA from a GTF annotation.

    Create a new GTF with the earliest start and latest finish associated with each gene_id.

    Args:
        gtfs: A list of GTF pandas.DataFrames, imported using gtf2pandas() from this library.
        filter(bool): Remove pre-mRNA models for single-model single-exon genes.
                    This reduces model inflation/duplication.
    Returns:
        [pandas.DataFrame]: List of dataframes, one per file.
    """
    result = []
    for g in gtfs:
        # Aggregate.
        grped = g.groupby("parent_id")
        gnum = len(grped)
        pres = pd.DataFrame( data= {
            "chr" : grped.head(1)["chr"],
            "source" : ['based_on_Araport11'] * gnum,
            "feature" : ['exon'] * gnum,  # for the benefit of existing 3rd-party parsers
            "start" : grped["start"].min(),
            "stop" : grped["stop"].max(),
            "score" : ['.'] * gnum,
            "strand" : grped.head(1)["strand"],
            "phase" : ['.'] * gnum,
            "attributes" : grped.head(1)["parent_id"],  # Temporary value.
            "parent_id" : grped.head(1)["parent_id"]
        })
        # Filter
        if (filter):
            u = grped["target_id"].apply(lambda x: len(set(x)))  # Number of models per gene.
            e = grped["feature"].apply(lambda x: Counter(x)["exon"])  # Number of exons per gene.
            idx = pres.index
            for gid in idx:
                # Drop pre-mRNA entries for single-model single-exon genes.
                if (u[gid] == 1 and e[gid] == 1):
                    pres.drop(gid, axis=0, inplace=True)
        # Format the attributes.
        for gid in pres.index :
            pres.ix[gid, "attributes"] = 'transcript_id \"' + gid + '_pre\"; gene_id \"' + gid + '\";'
        # Order columns to match GTF specification.
        result.append(pres[["chr","source", "feature", "start", "stop", "score", "strand", "phase", "attributes"]])
    return result


def samPatternStats(pattern, bam='-', bco=-4, bcl=4, literal=True, mmCap=2, wild='N', minFreq=0.01, filtered=False):
    """Find a pattern's positions and flanking sequences in an uncompressed header-less SAM.

    Meant to be used to verify the position of a known spacer sequence and identify
    the demultiplexing barcodes adjacent to it.

    Args:
        pattern(str): A sequence literal or regex pattern.
        sam(str): A BAM file.
        bco(int): How many positions before the tracer's start (-n) or after the
                    tracer's end (+n) is the barcode (-4)?
        bcl(int): How many nucleotides long is the barcode (4)?
        literal(bool): Is the pattern a literal sequence (True)? If so, mismatch patterns will also be created.
        mmCap(int): Set upper limit to number of mismatch positions that are allowed (2).
        wild(char): What singular character(s) is/are used for unknown values (N).
        minFreq(int): Discard results occuring in fewer than 100 reads.
        filtered(bool): The return Counters are filtered to remove rare events.
    Returns:
        List of Counters:
                    [0] Total number of reads (int)
                    [1] Number of reads that matched the anchor (int)
                    [2] collections.Counter object for read lengths
                    [3] collections.Counter object with read count of tracer matches
                    [4] collections.Counter object with read count of barcodes found
                        The barcode sequence is 'out_of_range' when the barcode is hanging over either end of the read
                        due to a very shifted tracer position.
                    [5] collection.Counter object with wildcard count downstream of the adapter region.
    """
    patlen = len(pattern)   # if literal
    Lengths = Counter()
    reads = 0
    matched = 0
    Positions = Counter()
    Barcodes = Counter()
    Wilds = Counter()
    samin = pysam.AlignmentFile(bam, 'rb', check_sq=False)  # Checking for reference sequences in the header has to be disabled for unaligned BAM.
    p = None
    if not literal:
        p = re.compile(pattern)
    # Search the pattern line by line.
    for line in samin:
        reads = reads + 1
        seq = line.get_forward_sequence()
        seqlen = len(seq)
        Lengths.update( ["\t".join(['Length', str(seqlen), '.', '.'])] )
        whichMinDist = None
        hit = None
        if literal:
            # Crawl along the sequence. When there is an exact match, this should break out after a few iterations. If there is any mismatch, it will have to go to the end.
            minDist = patlen      # Start with the largest possible hamming distance of all-mismatches.
            for i in range(0, seqlen - patlen):
                dist = lev.hamming(pattern, seq[i:(i+patlen)])
                if dist < mmCap and dist < minDist:
                    minDist = dist
                    whichMinDist = i
                    if dist == 0:
                        break   # Can't get any better.
        else:
            m = pattern.search(seq)
            if m :
                whichMinDist = m.start()       # No tie-breaker if anchor matches more than once. Just take the first. Provide more explicit patterns to reduce this occurence.
                hit = m.group(0)
        if whichMinDist is not None:        # if there is a match
            patend = whichMinDist + patlen
            if patend - 1 <= seqlen:    # Make sure nothing hangs off the end
                matched = matched + 1
                if literal:
                    hit = seq[whichMinDist:patend]
                Positions.update( ["\t".join(['Anchor', str(minDist), hit, str(whichMinDist + 1)])] )
                    # ie. Anchor    2   TTCCAGCATNGCTCTNAAAC    11
                # Identify the barcode
                if bco is not None and bcl is not None:
                    pos = patend - 1 + bco if bco > 0 else whichMinDist + bco
                    bcend = pos + bcl
                    if pos >= 0 and (bcend - 1) <= seqlen:  # Make sure nothing hangs off the end
                        Barcodes.update( ["\t".join(['Barcode', '', seq[pos:bcend], str(pos + 1)])] )
                            # ie. Barcode   .    ACGT 7
                # Count wildcards downstream.
                guidePos = max(bcend if bcend else 0, patend)      # whichever comes last, the barcode or the spacer.
                waggr = 0
                for w in wild:      # allow more than one wildcard characters
                    waggr = waggr + seq.count(w, guidePos)
                Wilds.update( ["\t".join(['Wildcards', str(waggr), '', ''])] )
    samin.close()
    # Filter out rare events to keep output uncluttered.
    if filtered:
        for k in list(Lengths.keys()):      # List gets all the values of the iterator before I edit the dict. That way the iterator doesn't crash.
            if Lengths[k] / reads * 100 < minFreq:
                del Lengths[k]
        for k in list(Positions.keys()):
            if Positions[k] / reads * 100 < minFreq:
                del Positions[k]
        for k in list(Barcodes.keys()):
            if Barcodes[k] / reads * 100 < minFreq:
                del Barcodes[k]
        for k in list(Wilds.keys()):
            if Wilds[k] / reads * 100 < minFreq:
                del Wilds[k]
    return [reads, matched, Lengths.most_common(), Positions.most_common(), Barcodes.most_common(), Wilds.most_common()]


def main(args):
    # Organize arguments and usage help:
    parser = argparse.ArgumentParser(description="Utility tasks relevant to sequencing.\
     Be sure to read the pydoc as well, to get more details about each taks. For the most \
     part each runtime task is associated to a library function.")

    # Input/Output.
    parser.add_argument('INPUTTYPE', type=str, choices=['L','T','D','P'],
                                help=" Specify the type of the TARGETs: \
                                'T' = The actual input filess. \
                                'L' = Text file(s) listing the input files. \
                                'P' = Get list of input files from STDIN pipe. \
                                'D' = Input data directly from STDIN pipe. \
                                ('D' is compatible with only some of the functions)")
    parser.add_argument('TARGET', type=str, nargs='*',
                                help=" The targets, space- or comma-separated. Usually files. \
                                Look into the specific task details below for special uses. \
                                Do not specify with INPUTTYPE 'P' or 'D'.")
    parser.add_argument('-O','--out', type=str, nargs=3,
                                help=" Send individual outputs to individual files instead of \
                                merging them to STDOUT. Output files will be like \
                                <out[0]>/<out[1]>target<out[2]>")
    # Parameters.
    parser.add_argument('-L','--log', action='store_true',
                                help=" Log this command to ./commands.log.")
    parser.add_argument('-c','--comments', action='store_true',
                                help=" Include commented info to STDOUT or files. (Default don't include)")
    parser.add_argument('-C','--STDERRcomments', action="store_false",
                                help=" Do NOT show info in STDERR. (Default show)")
    parser.add_argument('-v','--verbose', action="store_true",
                                help="Include more details/fields/etc in the output. Task-dependent.")
    # Tasks.
    parser.add_argument('--StarFinalLogs', type=str, choices=['tab','wiki'],
                                help="Combine multiple final summary mapping logs by STAR \
                                into a single text table of the specified format. Compact or verbose.")
    parser.add_argument('--premRNA', type=str, choices=['a', 'f'],
                                help="Infer pre-mRNA coordinates from a GTF annotation. Returns a GTF \
                                of pre-mRNA transcripts, comprising of the earliest start and latest finish \
                                coordinates for each gene. ** Choice 'a' returns a pre-mRNA for every gene, whereas \
                                choice 'f' filters out genes with only one transcript model comprising of a single exon \
                                ** Compatible with 'D' as INPUTTYPE.")
    parser.add_argument('--t2g', type=str, choices=['header', 'nohead'],
                                help="Extract transcript-gene ID pairs from a GTF file. The value determines\
                                whether to print a column header line or not.")
    parser.add_argument('--samFltrRegs', type=str,
                                help="Filter a headerless SAM file stream according to a SAM header file \
                                that contains the desired group of regions. This works only with the 'D' INPUTTYPE. \
                                The input stream is typically the output of `samtools view <somefile.bam>`. \
                                The output is streamed to STDOUT, typically to be piped back to `samtools view -b`. \
                                The header file will be prepended to the output stream.")
    parser.add_argument('--samPatternStats', type=str, nargs=5,
                                help="Number and location of matches of the pattern in the reads of a SAM stream. \
                                Arguments: [1] (str) anchor sequence, [2] (int) mismatches allowed in the anchor (use 'None' if anchor is regex),\
                                [3] (char) wildcard character(s) (like 'N' for unknown nucleotides),\
                                [4] (int) barcode offset (+n downstream of match end, \\-n upstream of match start, \
                                escaping the minus sign is important), [5] (int) barode length.")
    params = parser.parse_args(args)


    # CALL DETAILS.
    if params.log:
        import mylogs
        mylogs.log_command()
#    if params.STDERRcomments:
#        sys.stderr.write(ml.paramstring())


    # INPUT.
    flist = None
    if params.INPUTTYPE == 'P':
        # Read files list from STDIN
        flist = fu.FilesList()
        for line in sys.stdin:
            name = line.rstrip("\n").split("\t")[0]
            if name != "":
                flist.append(name)
    elif params.INPUTTYPE == 'L':
        # Create the FilesList, by appending the contents of all provided lists.
        flist = fu.FilesList().populate_from_files(params.TARGET)
    elif params.INPUTTYPE == 'T':
        # Create the FilesList by supplying a direct list of files.
        flist = fu.FilesList(params.TARGET)
    elif params.INPUTTYPE == 'D':
        # Data will be read from STDIN. No files needed. Make an empty list.
        # Not all functions will switch to STDIN given this. Several will simply do nothing.
        flist = fu.FilesList()
    else:
        sys.exit(ml.errstring("Unknown INPUTTYPE."))

    # OUTPUT.
    outstream = sys.stdout
    outfiles = None
    outdir, outpref, outsuff = None, None, None
    if params.out:
        outdir = fu.expand_fpaths([params.out[0]])[0]
        outpref = params.out[1]
        outsuff = params.out[2]
        outfiles = fu.make_names(flist.aliases, (outdir, outpref, outsuff))

    ### TASKS ###


    # Combine STAR LOGS.
    if params.StarFinalLogs:
        # Do it.
        df = collect_starFinalLogs(flist, all=params.verbose)
        # Call details.
        if params.comments:
                sys.stdout.write(ml.paramstring())
        # Formatting choice.
        if params.StarFinalLogs == "wiki":
            table = "^ " + df.index.name + " ^ " + " ^ ".join(df.columns.values.tolist()) + " ^\n"
            for row in df.itertuples():
                table += "|  " + " |  ".join(row) + " |\n"
            sys.stdout.write(table)
        else:
            sys.stdout.write(df.to_csv(sep="\t", header=True, index=True))
        # Done.
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("collecting STAR final logs"))


    # Create PRE-MRNA GTF.
    elif params.premRNA:
        # Import data and calculate the result.
        gtfs = gtf2pandas(flist)
        result = gtf2premrna(gtfs, filter=(params.premRNA == 'f'))
        # I need the for loop to iterate at least once. Relevant for STDIN input, since I have no input files listed then.
        if flist == []:
            flist.append("<STDIN>")
        # Print the contents.
        for i, (myfile, myalias) in flist.enum():
            if outfiles:
                # Send to individual file instead of STDOUT.
                outstream = open(outfiles[i], 'w')
            try:
                result[i].to_csv(outstream, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
            except IOError:
                pass
            finally:
                if outfiles:
                    # Don't want to accidentally close STDOUT.
                    outstream.close()
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("creating pre-mRNA annotation"))


    # Extract transcript and gene ID PAIRS from GTF
    elif params.t2g:
        # Import GTF.
        gtfs = gtf2pandas(flist)
        # I need the for loop to iterate at least once. Relevant for STDIN input, since I have no input files listed then.
        if flist == []:
            flist.append("<STDIN>")
        # Print the contents.
        hdr=None
        if params.t2g == "header":
            hdr = True
        else:
            hdr = False
        for i, (myfile, myalias) in flist.enum():
            if outfiles:
                # Send to individual file instead of STDOUT.
                outstream = open(outfiles[i], 'w')
            gtfs[i].dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
            try:
                gtfs[i].iloc[:,9:11].drop_duplicates().sort_values(["parent_id","target_id"]).to_csv(outstream, sep='\t', header=hdr, index=False, quoting=csv.QUOTE_NONE)
            except IOError:
                pass
            finally:
                if outfiles:
                    # Don't want to accidentally close STDOUT.
                    outstream.close()
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("extracting gene/transcript ID pairs."))


    # FILTER a BAM file by REGION
    elif params.samFltrRegs:
        if not params.INPUTTYPE == 'D':
            sys.exit("The only allowed INPUTTYPE is 'D' for streaming of header-less SAM content.")
        # Get regions from SAM header file
        regf = open(params.samFltrRegs, 'r')
        regions = list()
        p = re.compile('\sSN:(\S+)')
        for line in regf:
            sys.stdout.write(line)
            m = p.search(line)
            if m:
                regions.append(m.group(1))
        regf.close()
        # Parse SAM stream and output only the matching lines.
        p = re.compile('.+?\t\S+\t(\S+)')
        for r in sys.stdin:
            m = p.match(r)
            if m and m.group(1) in regions:
                sys.stdout.write(r)
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("filtering regions in SAM stream"))


    # Adapter STATS from SAM stream
    # read length distribution, pattern position distribution, pattern match
    elif params.samPatternStats:
        # Do.
        for i,f in enumerate(flist):
            result = samPatternStats(pattern=params.samPatternStats[0], bam=f, mmCap=int(params.samPatternStats[1]),
                                    bco=int(params.samPatternStats[3]), bcl=int(params.samPatternStats[4]),
                                    literal=(params.samPatternStats[1] != "None"), wild=params.samPatternStats[2], filtered=False)
            if outfiles:
                # Send to individual file instead of STDOUT.
                outstream = open(outfiles[i], 'w')
            # Print. The output should be an almost-tidy tab-delimited table.
            outstream.write( "\t".join(["Reads", '', '', '', str(result[0]), '(100% total)', os.path.basename(f) + "\n"]) )
            for v,c in result[2]:
                outstream.write( "\t".join([v, str(c), '(' + "{:.2f}".format(c / result[0] * 100) + '% total)', os.path.basename(f) + "\n"]) )
            outstream.write( "\t".join(["Matched", '', '', '', str(result[1]), '(' + "{:.2f}".format(result[1] / result[0] * 100) + '% total)', os.path.basename(f) + "\n"]) )
            for v,c in result[3]:
                outstream.write( "\t".join([v, str(c), '(' + "{:.2f}".format(c / result[0] * 100) + '% total)', os.path.basename(f) + "\n"]) )
            for v,c in result[4]:
                outstream.write( "\t".join([v, str(c), '(' + "{:.2f}".format(c / result[0] * 100) + '% total)', os.path.basename(f) + "\n"]) )
            for v,c in result[5]:
                outstream.write( "\t".join([v, str(c), '(' + "{:.2f}".format(c / result[0] * 100) + '% total)', os.path.basename(f) + "\n"]) )
            if outfiles:
                # Don't want to accidentally close STDOUT.
                outstream.close()
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("pattern stats in " + f))
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("pattern stats in all BAMs"))


#     # All done.
#     if params.STDERRcomments:
#         sys.stderr.write(ml.donestring())




#####    E X E C U T I O N   #####


# Call main only if the module was executed directly.
if __name__ == "__main__":
    main(sys.argv[1:])

    sys.exit(0)

#EOF
