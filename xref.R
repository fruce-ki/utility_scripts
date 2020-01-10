#!/usr/bin/env Rscript
    
library(data.table)

# Merge all files with a certain pattern form a folder, with a single other table, individually

args <- commandArgs(trailingOnly = TRUE)
mydir <- args[1]
suffix <- args[2]
xref <- args[3]
didx <- as.integer(args[4])
xidx <- as.integer(args[5])

# mydir <- "/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/M9262_slamseq/process/quant/"
# suffix <- "all_counts.txt"
# xref <- "/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/M9262_slamseq/aux/GRCm38/xref.tsv"
# didx <- 1
# xidx <- 1

if (is.na(didx)) didx <- 1
if (is.na(xidx)) xidx <- 1

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
