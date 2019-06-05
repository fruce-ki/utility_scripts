#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
# file <- '/Volumes/groups/obenauf/Kimon_Froussios/felix/MHC1_screen/process/pre_stats/CDMBHANXX_7_20190523B_20190526.tracer.stats'
outdir <- args[2]
# outdir <- ./process/pre_stats

library(tidyverse)

dt <- read_tsv(file,
               col_names = FALSE, na = '.', quoted_na = FALSE)
names(dt) <- c('type', 'score', 'seq', 'pos', 'reads', 'pc', 'bam')

total <- dt[[1, 'reads']]

bams <- unique(dt$bam)

for (b in bams) {
  subdt <- dt %>% filter(bam == b)
  if (b == '-') {
    b <- file   # If BAM was piped from STDIN, name this output after the stats output.
  } else {
    b <- file.path(outdir, b) # Typically I want output in a different directory than the BAM data
  }

  # Show the top stuff, save me scrolling through the stats
  write_tsv(dt %>% filter(reads / total > 0.001), paste(b, 'short.tsv', sep='.')) # Things greater than 0.1% of the data.

  # Summarise anchors by position
  anchorP <- dt %>%
    filter(type == 'Anchor') %>%
    group_by(pos) %>%
    summarise('type' = 'AnchorByPos', 'collectiveReads' = sum(reads)) %>%
    arrange(desc(collectiveReads))
  anchorP$collectivePc <- anchorP$collectiveReads / total * 100

  write_tsv(anchorP, paste(b, 'anchorPos.tsv', sep='.'))

  # Summarise anchors by mismatches
  anchorMM <- dt %>%
    filter(type == 'Anchor') %>%
    group_by(score) %>%
    summarise('type' = 'AnchorByMM', 'collectiveReads' = sum(reads)) %>%
    arrange(desc(collectiveReads))
  anchorMM$collectivePc <- anchorMM$collectiveReads / total * 100

  write_tsv(anchorMM, paste(b, 'anchorMM.tsv', sep='.'))

  # Summarise barcodes by sequence
  bcS <- dt  %>%
    filter(type == 'Barcode') %>%
    group_by(seq) %>%
    summarise('type' = 'BarcodeBySeq', 'collectiveReads' = sum(reads)) %>%
    arrange(desc(collectiveReads))
  bcS$collectivePc <- bcS$collectiveReads / total * 100

  write_tsv(bcS, paste(b, 'barcodeSeq.tsv', sep='.'))
}
