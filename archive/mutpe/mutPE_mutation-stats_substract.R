#!/usr/bin/env Rscript

# Substracts one (or more) stats files from another stats file.
# Individual mutation types at each position are summed, and the sums are then substracted across files.

# This script as-is WILL (probably) NOT WORK with the indel stats.

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("test", "/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/process_vdj/renamed/82001_CH12_C1r.aln_1-1.point.stats", "/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/process_vdj/renamed/82021_CH12_1_2.aln_1-1.point.stats")

outfile <- args[1]                            # Ouput stats file
fromfile <- args[2]                           # Input TSV stats to substract from (seq \t pos \t type \t count \t depth)
statsfiles <- args[3:length(args)]            # Input TSV stats to be substracted

if (file.exists(fromfile)){

  cntfilter <- 30
  covfilter <- 100
  
  # Left-hand side table to substract from
  lhs <- read_tsv(fromfile, col_types = 'cncnn')
  if (all(names(lhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
    names(lhs) <- c('chr', 'pos', 'type', 'count', 'depth')
  stopifnot(all(names(lhs) == c('chr', 'pos', 'type', 'count', 'depth')))
  
  # Only interested in collective mutation rate per position, not individual types.
  lhs <- lhs %>% 
    filter(grepl('^[ACGT]$', type)) %>%             # Getting reference bases, instead of mutated, ensures I get all the positions.
    mutate(type = "any", count = depth - count) %>%    # Aggregated mutation counts are the complement of reference counts.
    arrange(chr, pos) # %>% select(chr, pos, type, freq, depth)
  stopifnot(all(!is.na(lhs$pos)))
  
  
  # Substract cumulatively each right-hand side table
  # s <- statsfiles[1]
  for (s in statsfiles) {
    if (!file.exists(s)) {
      message(paste("No substraction file", s, ". Outputting", fromfile, "unchanged."))
      write_tsv(lhs, outfile)
      next
    }
    # Right-hand side table to be substracted in this round
    rhs <- read_tsv(s)
    if (all(names(rhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
      names(rhs) <- c('chr', 'pos', 'type', 'count', 'depth')
    stopifnot(all(names(rhs) == c('chr', 'pos', 'type', 'count', 'depth')))
    rhs <- rhs %>%
      filter(grepl('^[ACGT]$', type)) %>%             # Getting reference bases, instead of mutated, ensures I get all the positions.
      mutate(type = "any", count = depth - count) %>%    # Aggregated mutation counts are the complement of reference counts.
      arrange(chr, pos) %>%
      select(chr, pos, count, depth)
    stopifnot(all(!is.na(rhs$pos)))
  
    # Outter join, to fill in any potential gaps.
    lhs <- lhs %>%
      full_join(rhs, by=c('chr','pos'), fill)
    # Fix NAs from the merge, if any
    lhs$count.x[is.na(lhs$count.x)] <- 0
    lhs$count.y[is.na(lhs$count.y)] <- 0
    lhs$depth.x[is.na(lhs$depth.x)] <- 0
    lhs$depth.y[is.na(lhs$depth.y)] <- 0
    lhs$type[is.na(lhs$type)] <- "any"
    
    # Scale PER POSITION to lowest coverage and substract
    lhs <- lhs %>% 
      mutate(depth = pmin(depth.x, depth.y)) %>%
      # If the position is not covered at all in a track, I don't want to scale the other down to 0. Set to 1 instead.
      mutate(depth = pmax(1, depth)) %>%
      # If depth.x/.y is 0, then count.x/.y must be 0 too, so scaling factor is irrelevant, but still need to avoid division by 0.
      mutate(sf.x = ifelse(depth.x==0, 1, depth / depth.x), sf.y = ifelse(depth.y==0, 1, depth / depth.y)) %>%
      mutate(count.xs = count.x * sf.x, count.ys = count.y * sf.y) %>%
      mutate(count = count.xs - count.ys) %>%
      select(chr, pos, type, count, depth)
    stopifnot(all(!is.na(lhs$pos)))   # Something went wrong with the merge
  }
  
  write_tsv(lhs, outfile)

} else {
  message(paste("No base file", fromfile))
}





