#!/usr/bin/env Rscript

args         <- commandArgs(trailingOnly = TRUE)
file <- args[1]
# file <- '/Volumes/groups/obenauf/Kimon_Froussios/felix/MHC1_screen/process/pre_stats/CDMBHANXX_7_20190523B_20190526.tracer.stats'

library(tidyverse)

dt <- read_tsv(file,
               col_names = FALSE, na = '.', quoted_na = FALSE)
names(dt) <- c('type', 'score', 'seq', 'pos', 'reads', 'pc')

total = dt[[1, 'reads']]

# Show the top stuff, save me scrolling through the stats
write_tsv(dt %>% filter(reads > 10000), paste(file, 'short.tsv', sep='.'))

# Summarise anchors by position
anchorP <- dt %>%
  filter(type == 'Anchor') %>%
  group_by(pos) %>%
  summarise('type' = 'AnchorByPos', 'collectiveReads' = sum(reads)) %>%
  arrange(desc(collectiveReads))
anchorP$collectivePc <- anchorP$collectiveReads / total * 100

write_tsv(anchorP, paste(file, 'anchorPos.tsv', sep='.'))

# Summarise anchors by mismatches
anchorMM <- dt %>%
  filter(type == 'Anchor') %>%
  group_by(score) %>%
  summarise('type' = 'AnchorByMM', 'collectiveReads' = sum(reads)) %>%
  arrange(desc(collectiveReads))
anchorMM$collectivePc <- anchorMM$collectiveReads / total * 100

write_tsv(anchorMM, paste(file, 'anchorMM.tsv', sep='.'))

# Summarise barcodes by sequence
bcS <- dt  %>%
  filter(type == 'Barcode') %>%
  group_by(seq) %>%
  summarise('type' = 'BarcodeBySeq', 'collectiveReads' = sum(reads)) %>%
  arrange(desc(collectiveReads))
bcS$collectivePc <- bcS$collectiveReads / total * 100

write_tsv(bcS, paste(file, 'barcodeSeq.tsv', sep='.'))

