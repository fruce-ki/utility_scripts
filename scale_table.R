#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec <- matrix(c(
  'tableFile'     , 'i', 1, "character", "Tab-separated text table.",
  'outFile  '     , 'o', 1, "character", "Output file.",
  'help'          , 'h', 0, "logical"  , "Help.",
  'factor'        , 'v', 2, "numeric"  , "Scale all columns to this value (otherwise to smallest colSum)"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)
# opt <- list(tableFile='/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/HLA1_staggered_v2/process/crispr-processed/counts/library/counts_mageck_ctrls-grouped.txt', factor=50000000)

DT <- fread(opt$tableFile)
sel <- vapply(DT, class, character(1)) %in% c('numeric', 'integer')

lib <- colSums(DT[, sel, with=FALSE])
target <- NULL
if ( is.null(opt$factor) ) { 
  target <- min(lib)
} else {
  target <- opt$factor
}

fwrite(cbind( DT[, !sel, with=FALSE], 
              as.data.table(lapply (DT[, sel, with=FALSE], function(x) { x / sum(x) * target })) ),
       file = opt$outFile, sep="\t")

