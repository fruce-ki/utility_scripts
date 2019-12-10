#!/usr/bin/env Rscript
    
library(data.table)

# Merge all files with a certain pattern form a folder, with a single other table, individually

args <- commandArgs(trailingOnly = TRUE)
mydir <- args[1]
xref <- args[2]
didx <- args[3]
xidx <- args[4]
suffix <- args[5]

if (is.na(didx)) didx <- 0
if (is.na(xidx)) xidx <- 0

FF <- dir(mydir, full.names = TRUE)
FF <- FF[grep(suffix, FF)]

X <- fread(xref, colClasses="character")
# setkey(X, Entrez)
# length(unique(X$Entrez)) == length(unique(X$Name))
# length(unique(X$Entrez)) == length(unique(X$Id))

for (f in FF){
  DD <- fread(f, colClasses="character")
  # setkey(DD, Name)

  EE <- merge(X, DD, by.x=names(X)[xidx], by.y=names(DD)[didx])
  fwrite(file=sub('.txt$|.tsv$', '.xref.txt', f), EE, sep="\t")
}
