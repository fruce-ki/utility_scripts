#!/usr/bin/env Rscript

# Sum same-named columns across multiple tables.

library(data.table)

args <- commandArgs(trailingOnly=TRUE)
#args <- c('id', '/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/HLA1_staggered_v2/process/crispr-processed/counts/library/counts_mageck.txt', '/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/HLA1_staggered_v2/process/crispr-processed/counts/library/___counts_mageck.txt', '/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/HLA1_staggered_v2/description/borrowed_hilo.txt')
keycol <- args[1]         # Name of the row key column (must be shared among tables)
outfile <- args[2]        # output file
basefile <- args[3]       # first input table
addfiles <- args[4:length(args)] # additional input table(s), to be added to the values of the first table.

basedf <- fread(basefile)

#addfile <- addfiles[1]
for (addfile in addfiles){
  adddf <- fread(addfile)  
  stopifnot(nrow(basedf) == nrow(adddf))
  stopifnot(all(basedf[[keycol]] == adddf[[keycol]]))
  
  commoncols <- setdiff(intersect(names(basedf), names(adddf)), keycol)
  print( commoncols )
  
  for (N in names(commoncols)[2:13]) {
    basedf[, N] <- basedf[, N, with=FALSE] + adddf[, N, with=FALSE]
  }
}

fwrite(basedf, file=outfile, sep="\t")
