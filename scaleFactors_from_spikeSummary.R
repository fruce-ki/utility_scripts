#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
# args <- c('/Volumes/groups/zuber/zubarchive/USERS/Kimon/anja/M10716_slamseq/process/slam/spike_summary_counts_relevant.txt')


# Presumes the following columns (from my Split_BAM_Spikein.py script): row_ID, ambiguous_counts, spike_counts, subject_counts
DT <- read.table(args[1], stringsAsFactors=FALSE, sep="\t", header=TRUE)


# Geometric Mean : https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

target <- gm_mean(DT$spike_counts)


DT["spike_factor"] <- target / DT$spike_counts  # Multiplication factor to apply so as to get all spikeins to the same target value.
DT["subject_factor"] <- 1 / DT["spike_factor"]  # Multiplication factor to apply to the subject library,

write.table(DT[, c("row_ID", "subject_factor")], file=sub('txt|csv|tsv', 'factors.txt', args[1], perl=TRUE), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
