library(GenomicRanges)
library(rtracklayer)
library(data.table)


args <- commandArgs(trailingOnly=TRUE)
# args <- c("/groups/busslinger/Kimon/genome_mm10/refGene.mm10.2018_1206.withIgtcr.final.annotation.IMP_IGTCR_only.gtf")

# Input
DT <- as.data.table(rtracklayer::import(args[1]))
DT <- DT[type=="exon", ]

# Single-exon transcripts are not applicable
DT[, nexon := .N, by="transcript_id"]
singleexon <- DT[nexon <= 1, ]
DT <- DT[nexon > 1, ]
message(paste(length(unique(singleexon[, transcript_id])), "transcripts have only one exon. No introns will be created for them."))

# Look for problematic transcripts annotated on multiple chromosomes. In case this gene absurdity extends to transcripts.
DT[, numchr := length(unique(seqnames)), by="transcript_id"]
multichr <- DT[numchr > 1, ]
DT <- DT[numchr == 1, ]
message(paste(length(unique(multichr[, gene_id])), "transcripts are annotated each on multiple chromosomes. No introns will be created for them."))


# The next operation is order sensitive.
setorder(DT, transcript_id, start)

# Use exon ends to define intron ends.
DT[, irn_start := shift(end, n=1, type="lag") + 1, by='transcript_id']
DT[, irn_end := start - 1]
# Get rid of orphan rows (there is 1 less intron than exons in each transcript)
DT <- DT[!is.na(irn_start), ]

# Update properties.
DT[, exon_number := as.character(1:.N), by="transcript_id"]
DT[, type := 'intron']
DT[, source := 'filled_in']
DT[, start := irn_start]
DT[, end := irn_end]
DT[, irn_start := NULL]
DT[, irn_end := NULL]
DT[, nexon := NULL]
DT[, numchr := NULL]
DT[, width := (end + 1) - start]  #  GTF range left and right closed
DT[, exon_id := paste0(transcript_id, '_irn', exon_number)]

# Sanity check
DT <- DT[width >= 1, ]  # It is invalid for the end to be less than the start.
# library(ggplot2)
# ggplot(DT, aes(x=width)) +
#   geom_histogram(bins=500)
message(paste(length(unique(DT[, gene_id])), "genes received intron annotations."))
message(paste(length(unique(DT[, transcript_id])), "transcripts received intron annotations."))

# Output
setorder(DT, seqnames, start, end)
rtracklayer::export(makeGRangesFromDataFrame(unique(DT), keep.extra.columns = TRUE),
       sub(".gtf", ".only_introns.gtf", args[1]))
