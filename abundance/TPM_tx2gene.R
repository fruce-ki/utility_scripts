#!/usr/bin/env Rscript

library(tximport)
library(GenomicFeatures)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/mm10_ncbi/GCF_000001635.26_GRCm38.p6_genomic.gff', '/Volumes/groups/zuber/zubarchive/USERS/Kimon/sarah/philip2017/process/quantify', '/Volumes/groups/zuber/zubarchive/USERS/Kimon/sarah/philip2017/process/count/gene_tpm.tsv', '/Volumes/groups/zuber/zubarchive/USERS/Kimon/sarah/philip2017/process/count/gene_tpm_scaled.tsv')

GFF <- args[1]   # Annotation GFF3.
dir <- args[2]   # Parent directory containing only kallisto sample outputs.
OUT <- args[3]   # Destination TPM file.
OUT2 <- args[4]  # Destination pseudocount (scaled TPM) file.


txdb <- makeTxDbFromGFF(GFF)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

dups <- tx2gene$TXNAME[duplicated(tx2gene$TXNAME)]
ambiguous <- tx2gene[(tx2gene$TXNAME %in% dups),]
print(tx2gene[(tx2gene$GENEID %in% ambiguous$GENEID),])

#samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples <- list.files(dir)
files <- file.path(dir, samples, "abundance.h5")
all(file.exists(files))

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
txi.sum <- summarizeToGene(txi.kallisto, tx2gene)
colnames(txi.sum$abundance) <- samples
colnames(txi.sum$counts) <- samples

write.table(txi.sum$abundance, file=OUT, sep="\t", row.names = TRUE, quote = FALSE)
write.table(txi.sum$counts, file=OUT2, sep="\t", row.names=TRUE, quote = FALSE)
