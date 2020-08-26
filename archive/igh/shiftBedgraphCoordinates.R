#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
bedgraph <- args[1]
offset <- as.integer(args[2])
out <- args[3]

library(data.table)

dt <- fread(bedgraph)

dt[, V2 := V2 + offset]
dt[, V3 := V3 + offset]

fwrite(dt, file=out, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
