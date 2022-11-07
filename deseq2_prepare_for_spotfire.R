#!/usr/bin/env Rscript

library(getopt)
library(data.table)
library(purrr)

spec = matrix(c(
  'de',            'D', 1, "character", "Full path to Deseq2 results (aggregated or individual).",
  'cnt',           'C', 1, "character", "Full path to featureCounts table (aggregated or not), for the gene positions.",
  'fcThresh',      'f', 1, "numeric", "Fold change threshold.",
  'pCutoff',       'p', 1, "numeric", "P-value cut-off.",
  'countThresh',   'c', 1, "integer", "Count threshold upper cap.",
  'tpmThresh',     't', 1, "numeric", "TPM threshold upper cap.",
  'simplify',      's', 0, "logical", "Simplify column names, by removing origin, subset and formula, as long as this does not create duplicate names.",
  'url',           'u', 1, "character", "Genome Browser session URL to which the position will be appended for each gene."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(de='/Volumes/groups/busslinger/Kimon/sarah/R13546_RNAseq/results/DE/intron_genecounts/intron_genecounts.type_ctrl.LRT_cell.deseq2.tsv', cnt='/Volumes/groups/busslinger/Kimon/sarah/R13546_RNAseq/process/featureCounts/intron_genecounts.txt', fcThresh=0.5, pCutoff=0.0001, countThresh=500L, tpmThresh=25, simplify=TRUE, url="foobar")

if (is.null(opt$de)) stop("No input specified.")
if (is.null(opt$cnt) & !is.null(opt$url)) stop("Need a counts table for gene positions with which to form the URLs.")
if (is.null(opt$fcThresh)) { opt$fcThresh <- 2 } else {  opt$fcThresh <- as.numeric(opt$fcThresh)  }
if (is.null(opt$pCutoff)) { opt$pCutoff <- 0.05 } else { opt$pCutoff <- as.numeric(opt$pCutoff) }
if (is.null(opt$countThresh)) { opt$countThresh <- 100 } else { opt$countThresh <- as.numeric(opt$countThresh) }
if (is.null(opt$tpmThresh)) { opt$tpmThresh <- 50 } else { opt$tpmThresh <- as.numeric(opt$tpmThresh) }
if (is.null(opt$simplify)) opt$simplify <- FALSE


### Counts

# only used to extract the gene coordinates

CNT <- fread(opt$cnt)[, 1:4]  # featureCounts

CNT[, schr := lapply(strsplit(Chr, ';', fixed=TRUE), function(x){which(!grepl('_', x))})]  # drop alternative contigs etc and keep only main assembly
CNT[, Chr := map2(strsplit(Chr, ';', fixed=TRUE), schr, function(a, b){a[b]})]
CNT[, Start := map2(strsplit(Start, ';', fixed=TRUE), schr, function(a, b){a[b]})]
CNT[, End := map2(strsplit(End, ';', fixed=TRUE), schr, function(a, b){a[b]})]
CNT[, altonly := vapply(Chr, length, integer(1))==0]                            # Genes only present in alternative contigs and not the main assembly.
CNT[, scattered := vapply(Chr, function(x){length(unique(x))>1}, logical(1))]   # Find genes with transcripts annotated to multiple main chromosmes
CNT[(!scattered) & !altonly, chr := vapply(Chr, unique, character(1))]
CNT[!is.na(chr), start := vapply(Start, function(x){min(as.integer(unlist(x)))}, integer(1))]
CNT[!is.na(chr), end := vapply(End, function(x){max(as.integer(unlist(x)))}, integer(1))]
CNT[!is.na(chr), URL := paste0(opt$url, '&position=', chr, '%3A', start, '-', end)]               # %3a == :
CNT[(scattered), URL := "Ambiguous annotation to multiple chromosomes."]
CNT[(altonly), URL := "Not annotated on the main contigs."]
CNT <- CNT[, .(Geneid, URL)]


### DE

DE <- fread(opt$de)


# Identify columns
if (grepl('_vs_|all', opt$de)) {
  # pairwise Wald test
  lfc <- names(DE)[which(grepl("log2FoldChange(?!.shrink)", perl=TRUE, names(DE)))]
  lfcs <- names(DE)[which(grepl("log2FoldChange.shrink", names(DE)))]
  fc <- NULL
  p <- names(DE)[which(grepl("padj", names(DE)))]
  pv <- names(DE)[which(grepl('pvalue', names(DE)))]
  mlp <- names(DE)[which(grepl("mlog10p", names(DE)))]
  cnt <- names(DE)[which(grepl("maxCount", names(DE)))]
  scnt <- names(DE)[which(grepl("maxScaledCount", names(DE)))]
} else if (grepl('\\.LRT_', opt$de)) {
  # Likelihood Ratio test, probably not be pairwise
  lfc <- NULL
  lfcs <- NULL
  fc <- names(DE)[which(grepl("FC$", names(DE)))]
  p <- names(DE)[which(grepl("padj", names(DE)))]
  pv <- names(DE)[which(grepl('pvalue', names(DE)))]
  mlp <- names(DE)[which(grepl("mlog10p", names(DE)))]
  cnt <- names(DE)[which(grepl("baseMean", names(DE)))]
  scnt <- names(DE)[which(grepl("baseScaled", names(DE)))]
} else {
  stop("Unable to intepret contents from filename.")
}

## Create filters
#################

# for (X in lfc) {
#   # X <- lfc[1]
#   newcol <- sub("log2FoldChange", "FC_thresh", X)
#   set(DE, i=NULL, j=newcol, value=round(2 ^ abs(DE[[X]]), 1) )  #  cancel direction and lose excess decimals
# }

# FC threshold
for (X in lfcs) {
  # X <- lfcs[1]
  newcol <- sub("log2FoldChange.shrink", "sFC_thresh", X)
  set(DE, i=NULL, j=newcol, value="small")     # default value
  steps <- unique(c(2, 3, opt$fcThresh))
  steps <- steps[order(steps)]                      # smaller to bigger, order is important
  for (Y in steps)                                  # overwrite
    set(DE, i=which(2 ^ abs(DE[[X]]) >= Y), j=newcol, value=paste(">=", Y))
}

for (X in fc) {
  # X <- fc[1]
  newcol <- sub("FC", "FC_thresh", X)
  set(DE, i=NULL, j=newcol, value="small")     # default value
  steps <- unique(c(2, 3, opt$fcThresh))
  steps <- steps[order(steps)]                      # smaller to bigger, order is important
  for (Y in steps)                                  # overwrite
    set(DE, i=which(DE[[X]] >= Y | DE[[X]] < (1 / Y)), j=newcol, value=paste(">=", Y))
}


# Significance cutoff
for (X in p) {
  # X <- p[1]
  newcol <- sub("padj", "p_cutoff", X)
  set(DE, i=NULL, j=newcol, value="non-sig.")     # default value
  steps <- unique(c(0.05, opt$pCutoff))
  steps <- steps[order(steps, decreasing=TRUE)]     # high to low, order is important
  for (Y in steps)                                  # overwrite
    set(DE, i=which(abs(DE[[X]]) < Y), j=newcol, value=paste("<", Y))
}

# Global correction
P <- unique(DE[, c('name', pv), with=FALSE])
P <- melt(P, id.vars='name', value.name='pval', variable.name='headline')
P[, gpadj := p.adjust(pval, 'BH')]
P[, headline := sub('pvalue', 'global_padj', headline)]
P <- dcast(P, name ~ headline, value.var='gpadj')
DE <- merge(DE, P, by='name', all.x=TRUE)

# Global significance, cutoff
p <- names(DE)[which(grepl("global_padj", names(DE)))]
for (X in p) {
  # X <- p[1]
  newcol <- sub("padj", "p_cutoff", X)
  set(DE, i=NULL, j=newcol, value="non-sig.")     # default value
  steps <- unique(c(0.05, opt$pCutoff))
  steps <- steps[order(steps, decreasing=TRUE)]     # high to low, order is important
  for (Y in steps)                                  # overwrite
    set(DE, i=which(abs(DE[[X]]) < Y), j=newcol, value=paste("<", Y))
}

# FC direction
for (X in lfc) {
  # X <- lfc[1]
  newcol <- sub("log2FoldChange", "Deregulation", X)
  padj <- sub("log2FoldChange", "global_padj", X)
  set(DE, i=NULL, j=newcol, value="neutral")
  set(DE, i=which(DE[[X]] > 0 & DE[[padj]] < opt$pCutoff), j=newcol, value="UP")
  set(DE, i=which(DE[[X]] < 0 & DE[[padj]] < opt$pCutoff), j=newcol, value="DOWN")
}

# Count threshold cap
for (X in cnt) {
  # X <- cnt[1]
  newcol <- sub('maxCount|baseMean', 'count_thresh', X)
  set(DE, i=NULL, j=newcol, value=DE[[X]])
  set(DE, i=which(DE[[X]] > opt$countThresh), j=newcol, value=opt$countThresh) # Set upper cap
  set(DE, i=NULL, j=newcol, value=floor(DE[[X]]) )                             # Decimal precision is not useful for this
  if (grepl('maxCount', X))
    set(DE, j=X, value=NULL)
}

# TPM threshold cap
for (X in scnt) {
  # X <- scnt[1]
  newcol <- sub('maxScaledCount|baseScaled', 'scaledCount_thresh', X)
  set(DE, i=NULL, j=newcol, value=DE[[X]])
  set(DE, i=which(DE[[X]] > opt$tpmThresh), j=X, value=opt$tpmThresh) # Set upper cap
  set(DE, i=NULL, j=X, value=round(DE[[X]], 1) )                      # One decimal should be enough. Integers may be already too big for deep libraries.
}

## Spotfire does no recognise Inf values. So replace with suitable finite values.
# This affects -log(p) because DESeq2 can assign 0 to p.
for (X in mlp) {
  # X <- mlp[1]
  set(DE, i=which(is.infinite(DE[[X]])), j=X, value=max(60, max( DE[[X]][is.finite(DE[[X]])] )) )
}
for (X in fc) {
  # X <- fc[1]
  if (grepl('_vs_', opt$de)) {
    set(DE, i=which(is.infinite(DE[[X]]) & DE[[X]] > 0), j=X, value=max(10, max( DE[[X]][is.finite(DE[[X]])] )) )
    set(DE, i=which(is.infinite(DE[[X]]) & DE[[X]] < 0), j=X, value=min(10, min( DE[[X]][is.finite(DE[[X]])] )) )
  } else {
    set(DE, i=which(is.infinite(DE[[X]]) & DE[[X]] > 1), j=X, value=max(1000, max( DE[[X]][is.finite(DE[[X]])] )) )
    set(DE, i=which(DE[[X]] == 0), j=X, value=min(0.001, min( DE[[X]][DE[[X]] > 0] )) )
  }
}


# Inf fold-change is already handled by DESeq2 by default, by extrapolating a more likely FC from a fitted model and by explicitly capping the FC.


# Simplify column names
if (opt$simplify) {
  newnames <- sub("[^.]*?counts\\.", "", names(DE))   # counting mode
  if (all(!(duplicated(newnames))))
    setnames(DE, newnames)
  newnames <- sub("~[^.]*?\\.", "", names(DE))        # DE formula
  if (all(!(duplicated(newnames))))
    setnames(DE, newnames)
  newnames <- sub("\\.[^.]+$", "", names(DE))            # Contrast
  if (all(!(duplicated(newnames))))
    setnames(DE, newnames)
}

# Add URL
DE <- merge(DE, CNT, by.x='name', by.y='Geneid', all.x=TRUE, all.y=FALSE)


fwrite(DE, file=sub("txt$|tsv$", "spotfire.txt", opt$de), sep="\t", quote=FALSE, col.names=TRUE)


# ## Long form count maybe useful for some plots
# 
# CNT <- melt(CNT, id.vars="Geneid", measure.vars=names(CNT)[7:length(CNT)], variable.name="Sample", value.name="Count")
# TPM <- melt(TPM, id.vars="id", measure.vars=names(TPM)[2:length(TPM)], variable.name="Sample", value.name="TPM")
# TPM[, Sample := sub("\\..PM$", "", Sample)]
# setnames(TPM, c(names(CNT)[1:2], "TPM"))
# 
# CNT <- merge(CNT, TPM, all=TRUE)
# 
# CNT[, log10Count := log10(Count + 0.1)]
# CNT[, log10TPM := log10(TPM + (min(TPM[TPM>0], na.rm=TRUE)/10))]
# 
# 
# fwrite(CNT, 
# 			 file=sub("txt$", "spotfire.txt", opt$cnt), 
# 			 sep="\t", quote=FALSE, col.names=TRUE)
# 
