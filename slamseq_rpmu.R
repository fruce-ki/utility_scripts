#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/users/kimon.froussios/zuber/markus/OTI_slamseq/process/slam/dunk/count/SC1-A_lacxbeling-ctrl_4sU_no-ch_trimmed.fq_slamdunk_mapped_filtered_tcount.tsv')


for (f in args) {
  # f <- args[1]
  DT <- fread(f)
  
  # Add expression representations
  DT[, nonTcReadCount := ReadCount - TcReadCount]
  DT[, RPM := ReadCount / sum(ReadCount) * 1e6]
  DT[, RPMu := TcReadCount / sum(nonTcReadCount) * 1e6]
  
  # A gene can have multiple UTR features annotated. All quantified separately. So Entrez is not a unique ID.
  DT[, rowID := paste(Name, Chromosome, Start, End, Strand, sep="_")]
  
  # Reorder, and dump the extras.
  DT <- DT[, .(rowID, Name, ConversionRate, ReadsCPM, Tcontent, CoverageOnTs, ConversionsOnTs, ReadCount, TcReadCount, multimapCount, ConversionRateLower, ConversionRateUpper, RPM, RPMu)]
  # Prepare column headers for merging
  newnames <- paste( sub('_trimmed.fq_slamdunk_mapped_filtered_tcount.tsv', '', basename(f), fixed=TRUE), names(DT), sep="." )
  newnames[c(1,2)] <- names(DT)[c(1,2)]
  setnames(DT, newnames)
  
  fwrite(DT, file=sub('\\.[A-Za-z]+$', '.rpmu.txt', f, perl=TRUE), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=FALSE)
}

