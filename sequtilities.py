#!/sw/opt/python/2.7.3/bin/python

"""sequtilities.py

Author: Kimon Froussios
Date last revised: 20/09/2016
Python tested: 2.7.3

Library of utility functions relevant to sequencing tasks.
"""


import sys, argparse, re, csv
import pandas as pd

import fileutilities as fu
import mylogs as ml


def collect_starFinalLogs(flist, all=False):
    """Combine the listed Log.final.out files into a pandas Dataframe.
    
    File identifiers (filenames) will be trimmed at '_Log'.
    
    Args:
        flist: A FilesList or list of input files.
        all(bool): Show all fields (False).
    Returns:
        pandas.DataFrame
    """
    rows = None
    if all:
        # Still discording some irrelevant stuff and blank columns.
        rows = [2,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,27,28,29]
    else:
        rows = [4, 5, 8, 9, 23, 25, 27, 28, 29]
    df = fu.get_crosspoints(flist, cols=[1], rows=rows, colSep=["|"], header=False, index=0, merge=True)[0]
    
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


def gtf2premrna(flist):
    """Infer pre-mRNA from a GTF annotation.
    
    Args:
        flist: A FilesList or list of input GTF files.
    Returns:
        [pandas.DataFrame]: List of dataframes, one per file.
    """
    # Use my own parser, because it already handles input from multiple files or STDIN.
    input = fu.get_columns(flist, cols=range(0,9), colSep=["\t"], header=False, index=None, merge=False)
    result = []
    for i, gtf in enumerate(input):
        gtf.columns = ["chrom", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"]
        gtf["gene"] = gtf["attributes"].str.extract('gene_id \"?(\w+)')
        gtf.set_index(gtf["gene"], inplace=True)
        # Aggregate.
        grped = gtf.groupby("gene")
        gnum = len(grped)
        pres = pd.DataFrame( data= {
                "chrom" : grped.head(1)["chrom"],
                "source" : ['based_on_Araport11'] * gnum,
                "feature" : ['pre-mRNA'] * gnum,
                "start" : grped["start"].min(),
                "stop" : grped["stop"].max(),
                "score" : ['.'] * gnum,
                "strand" : grped.head(1)["strand"],
                "phase" : ['.'] * gnum,
                "attributes" : grped.head(1)["gene"],  # Temporary value.
                "gene" : grped.head(1)["gene"]
        })
        # Format attributes.
        for j in range(0, gnum) :
            pres.ix[j, "attributes"] = 'transcript_id \"' + pres.ix[j, "gene"] + '_pre\"; gene_id \"' + pres.ix[j, "gene"] + '\";'
        # Order columns to match GTF specification.
        result.append(pres[["chrom","source", "feature", "start", "stop", "score", "strand", "phase", "attributes"]])
    return result


def main(args):
    # Organize arguments and usage help:
    parser = argparse.ArgumentParser(description="Utility tasks relevant to sequencing.\
     Be sure to read the pydoc as well, to get more details about each taks. For the most \
     part each runtime task is associated to a library function.")

    # Input/Output.
    parser.add_argument('INPUTTYPE', type=str, choices=['T','L','P','D'],
                                help="Specify the type of the TARGETs: \
                                'T' = The actual targets, \
                                'P' = Pipe the list of targets from STDIN, \
                                'L' = Text file(s) listing the targets \
                                'D' = Pipe data from STDIN instead of from files. \
                                ('D' is compatible with only some of the functions)")
    parser.add_argument('TARGET', type=str, nargs='*',
                                help="File(s). ** Mandatory unless INPUTTYPE is 'P'.")
    parser.add_argument('-O','--out', type=str, nargs=3,
                                help="Send individual outputs to individual files instead of \
                                merging them and outputting to STDOUT. The first value is the \
                                target directory, the second is a common prefix, and \
                                the third value is a common suffix/extension. Output files will be like \
                                <out[0]>/<out[1]>target<out[2]>, where 'target' is each value of TARGET. \
                                Prefix and suffix may be empty strings: \"\".")
    # Parameters.
    parser.add_argument('-L','--log', action='store_true',
                                help="Log this command.")
    parser.add_argument('-c','--comments', action='store_true',
                                help="Include commented info to STDOUT or files.")
    parser.add_argument('-C','--STDERRcomments', action="store_false",
                                help="Do NOT show progress info in STDERR.")
    parser.add_argument('-v','--verbose', action="store_true",
                                help="Include more details/fields/etc in the output. Task-dependent.")
    # Tasks.
    parser.add_argument('--StarFinalLogs', type=str, choices=['tab','wiki'],
                                help="Combine multiple final summary mapping logs by STAR \
                                into a single text table of the specified format. Compact or verbose.")
    parser.add_argument('--premRNA', action='store_true',
                                help="Infer pre-mRNA coordinates from a GTF annotation. Returns a GTF \
                                of pre-mRNA transcripts, comprising of the earliest start and latest finish \
                                coordinates for each gene. ** Compatible with 'D' as INPUTTYPE.")
    params = parser.parse_args(args)
    
    
    # CALL DETAILS.
    if params.log:
        import mylogs
        mylogs.log_command()
    if params.STDERRcomments:
        sys.stderr.write(ml.paramstring())
        
    
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
        flist = FilesList()
    else:
        sys.exit(ml.errstring("Unknown INPUTTYPE."))
    
    # OUTPUT.
    outdir, outpref, outsuff = None, None, None
    if params.out:
        outdir = expand_fpaths([params.out[0]])[0]
        outpref = params.out[1]
        outsuff = params.out[2]
        
        
    ### TASKS ###
    
    
    # Combine STAR LOGS.
    if params.StarFinalLogs:
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
            sys.stderr.write(ml.donestring("collect STAR final logs"))
    
    
    # Create PRE-MRNA GTF.
    elif params.premRNA:
        result = gtf2premrna(flist)
        # Create output filenames, if applicable. If [], then STDOUT.
        outfiles = fu.make_names(flist.aliases, (outdir, outpref, outsuff))
        outstream = sys.stdout
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
            try:
                sys.stderr.write(ml.donestring("create pre-mRNA annotation"))
            except IOError:
                pass
            



#     # All done.     
#     if params.STDERRcomments:
#         sys.stderr.write(ml.donestring())   
    



#####    E X E C U T I O N   #####


# Call main only if the module was executed directly.
if __name__ == "__main__":
    main(sys.argv[1:])
    
    sys.exit(0)
    
#EOF