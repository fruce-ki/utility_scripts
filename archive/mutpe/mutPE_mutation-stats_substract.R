#!/usr/bin/env Rscript

# Substracts one (or more) stats files from another stats file.
# Individual mutation types at each position are summed, and the sums are then substracted across files.

args <- commandArgs(trailingOnly = TRUE)
# args <- c("/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process/in-vivo/HDR2/86762-86754_C2_GCB1-Ge1_3-10.substracted.stats", "/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process/in-vivo/HDR2/86762_B18_C2_GCB1.aln_3-10.stats", "/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process/in-vivo/HDR2/86754_B18_C2_G_e1.aln_3-10.stats")

outfile <- args[1]                            # Ouput stats file
fromfile <- args[2]                           # Input TSV stats to substract from (seq \t pos \t type \t count \t depth)
statsfiles <- args[3:length(args)]            # Input TSV stats to be substracted

library(dplyr)
library(readr)

# Left-hand side table to substract from
lhs <- read_tsv(fromfile)
if (all(names(lhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
  names(lhs) <- c('seq', 'pos', 'type', 'count', 'depth')
stopifnot(all(names(lhs) == c('seq', 'pos', 'type', 'count', 'depth')))
# Only interested in collective mutation rate per position, not individual types.
lhs <- lhs %>%
  filter(grepl('^[ACGT]$|ins', type)) %>%             # Getting reference bases, instead of mutated, ensures I get all the positions. Insertions are at offset positions, so need to be caught explicitly.
  mutate(type = "any", freq = 1 - count / depth) %>%  # Mutation frequencies are the complement of reference frequencies.
  mutate(depth = 1) %>%                               # This allows re-using the visualisation script with pre-calculated frequencies instead of counts and depths.
  arrange(seq, pos)
stopifnot(all(!is.na(lhs$pos)))


# Substract cumulatively each right-hand side table
#s <- statsfiles[1]
for (s in statsfiles) {
  # Right-hand side table to be substracted in this round
  rhs <- read_tsv(s)
  if (all(names(rhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
    names(rhs) <- c('seq', 'pos', 'type', 'count', 'depth')
  stopifnot(all(names(rhs) == c('seq', 'pos', 'type', 'count', 'depth')))
  rhs <- rhs %>%
    filter(grepl('^[ACGT]$|ins', type)) %>%
    mutate(freq = 1 - count / depth) %>%
    arrange(seq, pos) %>%
    select(seq, pos, freq)
  stopifnot(all(!is.na(rhs$pos)))

  # Outer join, to deal with the possibility of missing positions or missing sequences in either table
  lhs <- lhs %>%
    full_join(rhs, by=c('seq','pos'))
  # Fix NAs from the merge, if any
  lhs$freq.x[is.na(lhs$freq.x)] <- 0
  lhs$freq.y[is.na(lhs$freq.y)] <- 0
  lhs$type[is.na(lhs$type)] <- "any"
  lhs$depth[is.na(lhs$depth)] <- "1"
  stopifnot(all(!is.na(lhs$pos)))   # Something went wrong with the merge
  
  # Substract
  lhs <- lhs %>%
    mutate(freq = freq.x - freq.y) %>%
    select(seq, pos, type, freq, depth)
}

names(lhs) <- c('seq', 'pos', 'type', 'count', 'depth')  # Fake the headers to comply with the expectations of the visualisation script.
write_tsv(lhs, outfile)
