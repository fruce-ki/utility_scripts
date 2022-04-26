#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
gtf <- args[1]
minlen <- args[2]
out <- args[3]

library(data.table)

dt <- fread(gtf)

len <- dt$V5 - dt$V4
sel <- len >= minlen

print(paste(sum(sel), "lines will be kept out of", dim(dt)[1]))

fwrite(dt[sel, 1:9], file=out, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
