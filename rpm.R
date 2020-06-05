#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec = matrix(c(
  'help'         , 'h', 0, "logical",   "Help",
  'countsFile'   , 'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples",
  'nidcols'      , 'i', 1, "numeric",   "Number of ID cols (at the start of the table"
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(countsFile='/groups/zuber/zubarchive/USERS/Kimon/markus/OTI_quantseq2/process/quant_dedup/counts/all_counts_xref.txt', nidcols=3)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if (is.null(opt$nidcols)) {
  opt$nidcols <- 1
}

DT <- fread(opt$countsFile)
totals <- colSums(DT[, (opt$nidcols + 1):length(DT)], na.rm=TRUE)
for( x in names(totals) ){
  # x <- names(totals)[1]
  DT[, c(x) := DT[[x]] / totals[x] * 1e6]
}

fwrite(DT, file=sub('.txt|.tsv', '_rpm.txt', opt$countsFile), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
