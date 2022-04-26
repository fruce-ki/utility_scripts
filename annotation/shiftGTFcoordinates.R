#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
gtf <- args[1]
offset <- as.integer(args[2])
out <- args[3]

library(data.table)

dt <- fread(gtf)

dt[, V4 := V4 + offset]
dt[, V5 := V5 + offset]

fwrite(dt[, 1:9], file=out, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
