#!/usr/bin/env Rscript

# Substracts one (or more) stats files from another stats file.
# Individual mutation types at each position are summed, and the sums are then substracted across files.

# This script as-is WILL (probably) NOT WORK with the indel stats.

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("test", "/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process_vdj/in-vivo/wt/86765_B18_GCB4.aln.point.stats", "/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process_vdj/in-vivo/wt/86759_B18_AIDKO_GCB2.aln.point.stats")

outfile <- args[1]                            # Ouput stats file
fromfile <- args[2]                           # Input TSV stats to substract from (seq \t pos \t type \t count \t depth)
statsfiles <- args[3:length(args)]            # Input TSV stats to be substracted

if (file.exists(fromfile)){
  # Left-hand side table to substract from
  lhs <- fread(fromfile)
  if (all(names(lhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
    setnames(lhs, c('chr', 'pos', 'type', 'count', 'depth'))
  stopifnot(all(names(lhs) == c('chr', 'pos', 'type', 'count', 'depth')))
  
  # Substract respective frequencies, one file at a time.
  for (s in statsfiles) {
    # s <- statsfiles[1]
    if (!file.exists(s)) {
      message(paste("No substraction file", s, "."))
      # write_tsv(lhs, outfile)
      next
    }
    
    # Mutation frequencies. 0 for reference positions. Keeping the reference for coverage info.
    lhs[, freq := count / depth]
    lhs[grepl('^[ACGTN]$', type, perl=TRUE), freq := 0]
    
    # Right-hand side table to be substracted in this round
    rhs <- fread(s)
    if (all(names(rhs) == c('amplicon', 'pos', 'type', 'freq', 'total')))
      names(rhs) <- c('chr', 'pos', 'type', 'count', 'depth')
    stopifnot(all(names(rhs) == c('chr', 'pos', 'type', 'count', 'depth')))
    # Mutation frequencies. 0 for reference positions. Keeping the reference for coverage info.
    rhs[, freq := count / depth]
    rhs[grepl('^[ACGTN]$', type, perl=TRUE), freq := 0]
    
    
    # Outter join, to fill in any potential gaps.
    lhs <- merge(lhs, rhs, by=c('chr','pos', 'type'), all=TRUE)
    # Fix NAs from the merge, if any
    lhs[is.na(freq.x), freq.x := 0]
    lhs[is.na(freq.y), freq.y := 0]
    lhs[is.na(depth.x), depth.x := 0]
    lhs[is.na(depth.y), depth.y := 0]
    # Keep most conservative depth
    lhs[, depth := min(depth.x, depth.y), by=pos]
    # Substract frequencies. Then scale new frequency up to the estimated count for the combined depth.
    # The substraction is only as good as the sample with the least amount of evidence.
    lhs[, count := (freq.x - freq.y) * depth ]
    # Cap negative pseudocounts at 0.
    lhs[count < 0, count := 0 ]
    
    # Lose the extra columns.
    lhs <- lhs[, .(chr, pos, type, count, depth)]
  }
  
  fwrite(lhs[, .(chr, pos, type, count, depth)], file=outfile, quote=FALSE, sep="\t")

} else {
  message(paste("No base file", fromfile))
}




