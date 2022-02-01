library(data.table)
library(GenomicRanges)
library(rtracklayer)


args <- commandArgs(trailingOnly=TRUE)
# args <- c("/groups/busslinger/Kimon/genome_mm10/refGene.mm10.2018_1206.withIgtcr.final.annotation.IMP_IGTCR_only.gtf")

# Input
DT <- as.data.table(rtracklayer::import(args[1]))
DT <- DT[type=="exon",]

# Single-exon transcripts are not applicable
DT[, nexon := length(unique(exon_id)), by=.(transcript_id, strand)] # yes, I have seen a gene annotated on both strands.
singleexon <- DT[nexon <= 1, ]
DT <- DT[nexon > 1, ]
message(paste(length(unique(singleexon[, gene_id])), "transcripts have only one exon."))

# Look for genes annotated on multiple chromosomes. They cause problems.
DT[, numchr := length(unique(seqnames)), by=gene_id]
multichr <- DT[numchr > 1, ]
DT <- DT[numchr == 1, ]
message(paste(length(unique(multichr[, gene_id])), "genes are annotated each on multiple chromosomes."))

# With those out of the way, min and max coordinate for each gene would normally be my naive proxy of an un-spliced RNA.
# DT[, gmin := min(start), by="gene_id"]
# DT[, gmax := max(end), by="gene_id"]
# DT[, glen := gmax - gmin]
# However, there are still genes that have exons/transcripts scattered across chromosomes that would result in massive pre-mRNAs of >1e7 nt length (mm10)!!
# library(ggplot2)
# ggplot(DT, aes(x = glen)) +
#   geom_histogram(bins=500) +
#   scale_y_log10()
# In normal genes, alternative transcripts typically have some overlap. In the above genes, that is often not the case.
DT[, ntrans := length(unique(transcript_id)), by=.(gene_id, strand)]
DT[, tmin := min(start), by=.(transcript_id, strand)]
DT[, tmax := max(end), by=.(transcript_id, strand)]
ovlpsearch <- unique(DT[, .(gene_id, ntrans, transcript_id, tmin, tmax, strand)])  # multi-exon transcripts cause false positives, so I gotta reduce it to one row per transcript.
setorder(ovlpsearch, gene_id, tmin, tmax, strand)
ovlpsearch[ntrans > 1, lagtmax := data.table::shift(tmax, n=1L), by=.(gene_id, strand)] #  want to compare consecutive transcripts
ovlpsearch[, overlap := lagtmax > tmin]
ovlpsearch[, numovlp := sum(overlap, na.rm=TRUE), by=.(gene_id, strand)]

message(paste(length(ovlpsearch[numovlp < ntrans - 1 & numovlp > 0, unique(gene_id)]), "genes have a mix of overlapping and non-overlapping transcripts. Assuming that a single unspliced RNA is still appropriate."))

scattered <- DT[gene_id %in% unique(ovlpsearch[numovlp == 0 & ntrans > 1, gene_id])]
DT <- DT[! gene_id %in% unique(scattered[, gene_id])]
message(paste(length(unique(scattered[, gene_id])), "genes have multiple transcripts scattered along a chromosome with no location overlaps at all among them."))

# The remaining genes should be safe to generate proxy pre-splice RNA by just min/max the gene's coordinates.
DT[, prestart := min(start), by=.(gene_id, strand)]
DT[, preend := max(end), by=.(gene_id, strand)]
DT[, prewidth := preend - prestart + 1]
DT[, presource := 'minmax']
DT[, pretranscript_id := paste0(gene_id, "_", strand, "pre")]
DT[, preexon_id := paste0(gene_id, "_", strand, "pre.1")]
DT <- DT[exon_number == 1, .(seqnames, prestart, preend, prewidth, strand, presource, type, score, phase, gene_id, pretranscript_id, exon_number, preexon_id, gene_name)]
setnames(DT, sub('pre', '', names(DT)))
message(paste(length(unique(DT[, gene_id])), "gene-wise unspliced RNA annotations created."))

# For the multi-chromosomal and scattered genes create per-transcript unspliced RNA
scattered <- rbind(scattered[, names(multichr), with=FALSE], multichr)
scattered[, nexon := .N, by="transcript_id"]
scattered <- scattered[nexon > 1, ]
scattered[, prestart := min(start), by="transcript_id"]
scattered[, preend := max(end), by="transcript_id"]
scattered[, prewidth := preend - prestart + 1]
scattered[, presource := 'minmax']
scattered[, pretranscript_id := paste0(transcript_id, "_pre")]
scattered[, preexon_id := paste0(transcript_id, "_pre.1")]
scattered <- scattered[exon_number == 1, .(seqnames, prestart, preend, prewidth, strand, presource, type, score, phase, gene_id, pretranscript_id, exon_number, preexon_id, gene_name)]
setnames(scattered, sub('pre', '', names(scattered)))
message(paste(length(unique(scattered[, gene_id])), "transcript-wise unspliced RNA annotations created for scattered genes."))

# Output
DT <- rbind(DT, scattered)
setorder(DT, seqnames, start, end)
rtracklayer::export(makeGRangesFromDataFrame(unique(DT), keep.extra.columns = TRUE),
       sub(".gtf", ".only_premRNA.gtf", args[1]))

