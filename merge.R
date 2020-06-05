#!/usr/bin/env Rscript

## Full/Outer merge of any number of tables. The first column is used as index.
## Output to STDOUT, be sure to redirect!

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('~/zuber/markus/OTI_quantseq2/aux/GRCm38/xref2.tsv', '~/zuber/markus/OTI_quantseq2/results/quant_dedup/DE/cell-bulk_time-at.~gene.gene_Kdm2b_vs_Control.deseq2.tsv', '~/zuber/markus/OTI_quantseq2/results/quant_dedup/DE/cell-bulk_time-at.~gene.gene_Kdm2b_vs_Rc3h1.deseq2.tsv', '~/zuber/markus/OTI_quantseq2/results/quant_dedup/DE/cell-bulk_time-at.~gene.gene_Rc3h1_vs_Control.deseq2.tsv')

stopifnot(length(args) > 1) 

# Parse all tables.
tables <- lapply(args, fread)

# Enforce that the first column shares name across tables.
n <- names(tables[[1]])[1]
for (tab in tables) {
  N <- names(tab)
  N[1] <- n
  setnames(tab, N)
  setkeyv(tab, n)
}

# Merge and dump to STDOUT
fwrite(Reduce(function (x, y) { merge(x, y, all=TRUE, by=n ) }, tables),
       file='',
       sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)