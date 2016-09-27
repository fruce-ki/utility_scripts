# utility_scripts

A collection of general-purpose code that makes my life easier.

## Contents

### fileutilities.py

Batch tasks of file-munging that I seem to need frequently enough to merit automation. Works as script or library.  (python 2.7.3)

Library and script:

1. Probe files for access/type/emptiness. (rather pointless, kept as record of heuristic for differentiating binary from text)
2. List & filter directories contents and assign aliases.
3. Soft-link listed files into destination.
4. Loop shell command over listed files.
5. Replace delimiting character (or any substring really) with new one.
6. Count columns in delimited files.
7. Extract columns from delimited files, even ones with irregular row lengths.
8. Extract random columns from delimited files (i.e. for bootstrapping).
9. Append columns to delimited files (merge delimited files horizontally).
10. Get the set of different values in any row or column of a delimited file (for fields with loads of observations but few values).

Library only:

1. The FilesList class.
2. Human-friendly string list sorting.
3. Get values at intersection of specific rows and columns of delimited texts.
4. Auto-suffix duplicates in a list of strings.
5. Auto-generate file names using a common target path, prefix and suffix and a list of strings.
5. First / Last lines of listed files. Obsolete, but kept as record of buffer access tricks.

A partial testing suite available in test_fileutilities.py.

### mylogs.py

Custom logging of commands and messages in consistent format. Works as library or command-wrapping script. (python 2.7.3)

Library:

1. Readable current time stamp.
2. Format various message strings (call parameters, end of task, info, warning, error).
3. Log call details of main script to logfile.
4. Log message to file.

Script:

1. Log and then execute any shell command.
2. Log a message.

### sequtilities.py

Custom utility functions related to Bioinformatics sequence analysis. Works as library or script. (python 2.7.3)

Library and script:

1. Aggregate and summarize the listed final logs of STAR.
2. Infer pre-mRNA coordinates from a GTF genome annotation.

Library only:

1. Import GTF into a pandas.DataFrame. Extract gene_id and transcript_id as new columns and set gene id as index.
