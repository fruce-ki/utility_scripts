#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
column <- as.integer(args[2])
offset <- as.numeric(args[3])

# file="/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/M9262_slamseq/aux/GRCm38/3UTR/mm10_refseq_ensembl_3UTR_fillup.gtf"
# column <- 4L
# offset <- 1L

DT <- fread(file)

set(DT, j=column, value= DT[[column]] + offset)

fwrite(DT, file=file, sep = "\t", quote = FALSE)
