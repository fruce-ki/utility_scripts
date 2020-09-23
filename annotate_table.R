#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly=TRUE)

# args <- c('/Volumes/groups/zuber/zubarchive/USERS/Kimon/UTR_forensics/aux/GCF_000001635.26_GRCm38.p6_genomic.bed', '1', '/Volumes/groups/zuber/zubarchive/USERS/Kimon/UTR_forensics/aux/GCF_000001635.26_GRCm38.p6_xref.txt', '1', '2')

# Annotate first table with specific field from second table using specified fields as cross-index.
# Output to STDOUT, be sure to redirect!

left <- fread(args[1], header=FALSE)
idxl <- as.integer(args[2])

idxr <- as.integer(args[4])
idxa <- as.integer(args[5])
right <- unique(fread(args[3], select=c(idxr, idxa), header=FALSE))
names(right) = c('K','V')
setkey(right, K)

left[, Annotated := right[left[[idxl]], V] ]

fwrite(left, file='', sep="\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
