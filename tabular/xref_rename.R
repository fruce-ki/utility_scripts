#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
columns <- args[1]  # Table to replace column names
rows <- args[2]     # Table to replave row names
tbl <- args[3]      # Table whose names to replace
out <- args[4]      # Destination file 




covars <- fread(columns)
names(covars)[1:2] <- c("old", "new")
setkey(covars, old)

xref <- fread(rows)
names(xref)[1] <- "key"
setkey(xref, key)

counts <- fread(tbl)
names(xref)[1] <- "key"
names(counts) <- vapply(names(counts), function(n){
                      columns$new[ which( vapply(columns$old, function(x){grepl(x, n, fixed=TRUE)}, logical(1)) ) ]
                  }, character(1))

names(counts)[2:length(counts)] <- covars[names(counts)[2:length(counts)], name]
fwrite(merge(counts, xref, by.x= "key", by.y="key"), file=otu, sep="\t")
