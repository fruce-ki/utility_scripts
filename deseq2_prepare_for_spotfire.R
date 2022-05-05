#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec = matrix(c(
  # 'cnt',           'C', 1, "character", "Full path to aggregated featureCount results.",
  'de',            'D', 1, "character", "Full path to aggregated Deseq2 results.",
  # 'tpm',           'T', 1, "character", "Full path to aggregated TPM.",
  'lfcThresh',     'f', 1, "numeric", "log2 Fold change threshold.",
  'pCutoff',       'p', 1, "numeric", "P-value cut-off.",
  # 'countThresh',   'c', 1, "integer", "Count threshold.",
  # 'tpmThresh',     't', 1, "numeric", "TPM threshold.",
  'simplify',      's', 0, "logical", "Simplify column names, by removing origin, subset and formula, as long as this does not create duplicate names."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(de='/scratch-cbe/users/kimon.froussios/tanja/R13065_RNAseq/results/DE/intron_genecounts/intron_genecounts.TF_Ikzf1.~Condition.Condition_Exp_vs_Ctrl.deseq2.tsv', lfcThresh=2, pCutoff=0.01, countThresh=50L, tpmThresh=5, simplify=TRUE)

if (is.null(opt$de))
  stop("No input specified.")
# if (is.null(opt$cnt) | is.null(opt$tpm))
# 	stop("No input specified.")
if (is.null(opt$lfcThresh))
  opt$fcThresh <- 1
if (is.null(opt$pCutoff))
  opt$pCutoff <- 0.05
# if (is.null(opt$countThresh))
#   opt$countThresh <- 50
# if (is.null(opt$tpmThresh))
#   opt$tpmThresh <- 5
if (is.null(opt$simplify))
  opt$simplify <- FALSE

### Counts

# CNT <- fread(opt$cnt)
# TPM <- fread(opt$tpm)

### DE

DE <- fread(opt$de)

## Calculate filtering vectors at some preset levels

# Identify columns
fc <- names(DE)[which(grepl("FC", perl=TRUE, names(DE)))]
lfcs <- names(DE)[which(grepl("log2FoldChange.shrink", names(DE)))]
p <- names(DE)[which(grepl("padj", names(DE)))]
mlp <- names(DE)[which(grepl("mlog10p", names(DE)))]

# Create multiple-choice style filters

for (X in fc) {
  # X <- fc[1]
  newcol <- sub("FC", "FC_thresh", X)
  set(DE, i=NULL, j=newcol, value=abs(DE[[X]]))
}

for (X in lfcs) {
  # X <- lfcs[1]
  newcol <- sub("log2FoldChange.shrink", "slFC_thresh", X)
  set(DE, i=NULL, j=newcol, value="unfiltered")     # default value
  steps <- as.character(unique(c(1, 2, 3, opt$lfcThresh)))
  steps <- steps[order(steps)]                      # smaller to bigger, order is important
  for (Y in steps)                                  # overwrite
    set(DE, i=which(abs(DE[[X]]) >= Y), j=newcol, value=paste(">=", Y))
}

for (X in p) {
  # X <- p[1]
  newcol <- sub("padj", "p_cutoff", X)
  set(DE, i=NULL, j=newcol, value="unfiltered")     # default value
  steps <- unique(c(0.05, opt$pCutoff))
  steps <- steps[order(steps, decreasing=TRUE)]     # high to low, order is important
  for (Y in steps)                                  # overwrite
    set(DE, i=which(abs(DE[[X]]) < Y), j=newcol, value=paste("<", Y))
}

## Spotfire does no recognise Inf values. So replace with something numeric.
# This affects -log(p) because DESeq2 can assign 0 to p.
for (X in mlp) {
  # X <- mlp[1]
  set(DE, i=which(is.infinite(DE[[X]])), j=X, 
      value=max(60, max( DE[[X]][is.finite(DE[[X]])] )) )
}
# Inf fold-change is already handled by DESeq2 by default, by extrapolating a more likely FC from a fitted model and by explicitly capping the FC.


## Simplify column names

if (opt$simplify) {
  newnames <- sub("[^.]*?counts\\.", "", names(DE))
  if (all(!(duplicated(newnames))))
    setnames(DE, newnames)
  newnames <- sub("~[^.]*?\\.", "", names(DE))
  if (all(!(duplicated(newnames))))
    setnames(DE, newnames)
  newnames <- sub("\\..*$", "", names(DE))
  if (all(!(duplicated(newnames))))
    setnames(DE, newnames)
}

## Identify fields that can be dropped, to de-clutter the lists of columns in spotfire

keep <- names(DE)[!grepl("lfcSE|pvalue|meetthresh|istop|category|StDev|meanCount|meanScaled", names(DE))]
DE <- DE[, keep, with= FALSE]

# DE <- merge(merge(DE, CNT[, c(1, 7:length(CNT)), with=FALSE], by.x='name', by.y='Geneid', all=TRUE), TPM, by.x='name', by.y='id', all=TRUE)


fwrite(DE, file=sub("txt$|nolab.tsv$", "spotfire.txt", opt$de), sep="\t", quote=FALSE, col.names=TRUE)


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
