#!/usr/bin/env Rscript

library(getopt)
library(DESeq2)
library(data.table)

spec = matrix(c(
  'bmF',           'B', 0, "logical",   "Disable baseMean filtering. (False)",
  'baseDir',       'b', 1, "character", "Base directory for everything (.)",
  'lfcthreshold',  'c', 1, "numeric",   "Log2 fold-change threshold (1)",
  'createID',      'C', 0, "logical",   "Non-summarized intron or exon counts by featureCounts. The first column is not unique IDs because genes are repeated.",
  'differential',  'd', 1, "character", "Name of factor along which to compare ('Condition').",
  'countsFile',    'f', 1, "character", "Tab-separated table of raw counts with `row_ID` and all the samples.",
  'specnorm',      'G', 1, "character", "Special normalisation. Gene ids that match this pattern will not be included for calculation of library sizes for TPM/RPM normalisation.",
  'idcol',         'i', 1, "numeric",   "ID column to use (1). The others will be removed. See also -I.",
  'nidcols',       'I', 1, "numeric",   "Number of ID columns at the start of the table (1). See also -i.",
  'prescaled',     'k', 0, "logical",   "Don't let Deseq2 rescale the libraries. (False)",
  'label',         'l', 0, "logical",   "Add comparison details to the standard column names. Useful for merging multiple outputs.",
  'minSingle',     'm', 1, "numeric",   "Minimum number of reads in any sinlge sample for a gene to be considered if the minMean is not met (100).",
  'minMean',       'M', 1, "numeric",   "Minimum mean number of reads across all samples combined for a gene to be considered (10).",
  'ntop',          'n', 1, "numeric",   "Number of hits to highlight (50)",
  'resultsDir',    'o', 1, "character", "Directory in which to save the contrast results. If omitted, only the RDS will be output.",
  'pcutoff',       'p', 1, "numeric",   "P-value cutoff (0.05)",
  'prefix',        'P', 1, "character", "Prefix to add to output file name and column names (''). ",
  'RDSoutdir',     'r', 1, "character", "Directory in which to save the raw DESeq2 object (./process)",
  'reducedFormula','R', 1, "character", "Reduced formula for LR test",
  'samplesFile',   's', 1, "character", "Tab-separated table with `Sample` column followed by the variable columns. Values must NOT contain punctuation other than '_' !!!",
  'reportTemplate','T', 1, "character", "Template Rmd file for DE report.",
  'vst',           'v', 0, "logical",   "Plot VST counts instead of rlog2",
  'widthsFile',    'w', 1, "character", "Full path to table of feature lengths",
  'widthsCol',     'W', 1, "integer",   "Column in the counts file that contains feature lengths. This should be among the non-count columns, see also -I. It is an alternative to providing a separate lengths file.",
  'comparisons',   'x', 1, "character", "Comma separated list indicating which levels of the formula factor to compare, in order of appearance. ie '2v1,3v1,3v2,4v5'."
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(createID=FALSE, baseDir='/scratch-cbe/users/kimon.froussios/sarah/R13546_RNAseq', countsFile='/process/featureCounts_fixed/exon_genecounts.txt', resultsDir='results/DE', RDSoutdir='process/DE', samplesFile='description/covars_de.txt', minMean=0, minSingle=0, lfcthreshold=1, nidcols=6, idcol=1, ntop=50, bmF=FALSE, pcutoff=0.05, prescaled=FALSE, prefix='', label=FALSE, widthsCol=6, reportTemplate="/groups/busslinger/Kimon/sarah/R13546_RNAseq/code/deseq2_report_template_SF.Rmd", differential='Condition', comparisons='2v1,4v3,6v5,8v7,10v9,12v11,14v13,16v15,18v17,20v19,22v21,24v23')

# opt <- list(createID=FALSE, baseDir='/Volumes/groups/busslinger/Kimon/sarah/R13546_RNAseq', countsFile='process/featureCounts/intron_genecounts.txt', resultsDir='results/DE', RDSoutdir='process/DE', samplesFile='description/covars_de.txt', minMean=0, minSingle=0, lfcthreshold=1, nidcols=6, idcol=1, ntop=50, bmF=FALSE, pcutoff=0.05, prescaled=FALSE, prefix='', label=FALSE, widthsCol=6, reportTemplate="/Volumes/groups/busslinger/Kimon/sarah/R13546_RNAseq/code/deseq2_report_template_SF.Rmd", differential='Condition', comparisons='2v1,4v3,6v5,8v7,10v9,12v11,14v13,16v15,18v17,20v19,22v21,24v23')


if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

stopifnot(!is.null(opt$countsFile))
stopifnot(!is.null(opt$samplesFile))

if (is.null(opt$reportTemplate))  opt$reportTemplate <- '~/utility_scripts/deseq2_report_template_SF.Rmd'
if (is.null(opt$createID))  opt$createID <- FALSE 
if ((!is.null(opt$specnorm)) && opt$specnorm == "NULL") opt$specnorm <- NULL
if (is.null(opt$prescaled)) opt$prescaled <- FALSE
if (is.null(opt$label)) opt$label <- FALSE
if (is.null(opt$baseDir)) opt$baseDir <- '.'
if (is.null(opt$RDSoutdir)) opt$RDSoutdir <- './process' 
if (is.null(opt$resultsDir))  opt$resultsDir <- './results' 
if (is.null(opt$minMean)) opt$minMean <- 10
if (is.null(opt$minSingle)) opt$minSingle <- 100
if (is.null(opt$ntop))  opt$ntop <- 50
if (is.null(opt$pcutoff)) opt$pcutoff <- 0.05
if (is.null(opt$lfcthreshold))  opt$lfcthreshold <- 1
if (is.null(opt$nidcols)) opt$nidcols <- 1
if (is.null(opt$idcol)) opt$idcol <- 1
if (is.null(opt$bmF)) opt$bmF <- FALSE
opt$bmF <- !opt$bmF    # Convert -B flag functionality from switch-on to switch-off. 
if (is.null(opt$vst)) opt$vst <- FALSE
if (is.null(opt$differential)) opt$differential <- 'Condition'


# Output destinations
dir.create(file.path(opt$baseDir, opt$RDSoutdir, opt$prefix), recursive=TRUE)
dir.create(file.path(opt$baseDir, opt$resultsDir, opt$prefix), recursive=TRUE)

# Covariates
covars <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", row.names="Sample", header=TRUE, check.names = FALSE, colClasses="character")

# Parse contrasts
conditions <- unique(covars[[opt$differential]])
comparisons <- strsplit(opt$comparisons, ',')[[1]]
comparisons <- lapply(comparisons, function(x){
  # x <- comparisons[1]
  x <- as.integer(strsplit(x, 'v')[[1]])
  c(opt$differential, conditions[x])
})

# Counts
cts <- fread(file.path(opt$baseDir, opt$countsFile), sep="\t", header=TRUE, check.names = FALSE, colClasses="character")
# Setting all columns to character prevents numeric IDs at the top of the ID columns from auto-defining these columns as numeric when there could be character IDs further down that would otherwise become NaN.
cts <- unique(cts)         # non-summarized featureCounts have repeated exons from sharing among transcripts.
xref <- cts[, 1:opt$nidcols]
cts <- as.matrix(cts[, rownames(covars), with=FALSE])    # Order and subset columns
cts <- matrix(as.numeric(cts), ncol=ncol(cts))
colnames(cts) <- row.names(covars)  # because the type conversion lost them
cts[is.na(cts)] <- 0

if (opt$createID) {
  xref[, id := paste0(Geneid, '.', Start, '_', End)]
  if (any(duplicated(xref$id)))
    stop()
  rownames(cts) <- xref$id
} else {
  rownames(cts) <- xref[[opt$idcol]]
}

# Feature lengths
if (!is.null(opt$widthsCol)) {
  featLens <- data.table(ID=rownames(cts),
                         LEN=as.integer(xref[[opt$widthsCol]]))
} else if (!is.null(opt$widthsFile)) {
  featLens <- fread(opt$widthsFile, sep="\t", header=TRUE, check.names = FALSE, colClasses=c("character", "numeric"))
  # Match order of feature lengths to order of features
  setkeyv(featLens, names(featLens)[1])
  featLens <- featLens[rownames(cts), ]
} else {
  featLens <- NULL
}
stopifnot(all(rownames(cts) == featLens[[1]]))


# Remove unused samples from the counts table,
# and ensure samples are in the same order as the covariates table.
cdn <- rownames(covars)
ctsn <- colnames(cts)
subcts <- cts[, ctsn[match(cdn, ctsn)]]
# Fill in non-detections
subcts[is.na(subcts)] <- 0

# DESeq2 begins here
DDS <- DESeqDataSetFromMatrix(countData = round(subcts, digits=0),
                              colData = covars,
                              design = as.formula(paste('~', opt$differential)))
# To scale or not to scale
if (opt$prescaled) {
  factors <- rep(1, times=length(cdn))
  names(factors) <- cdn
  sizeFactors(DDS) <- factors
}

# Prefilter out genes with very poor coverage. This will completely remove those genes from the output as well.
# DESeq2 internally applies additional thresholds. This only manifests as missing pvalues, the genes stay present.
keep <- rowMeans(counts(DDS)) >= opt$minMean | rowSums(counts(DDS) >= opt$minSingle) >= 1
DDS <- DDS[keep, ]

# Prepare file name.
if (opt$prefix == '') {
  prefix <- gsub('[ ~:]', '', opt$differential)
} else {
  prefix <- paste(opt$prefix, gsub('[ ~:]', '', opt$differential), sep='.')
}
ddsrds <- gsub('//', '/', file.path(opt$baseDir, opt$RDSoutdir, opt$prefix, paste0(prefix, '.deseq2data.RDS')))
name <- paste0(prefix, '.deseq2.html')

saveRDS(DDS, file=ddsrds)

#######################
#######################
# Fire up an Rmd report
rmarkdown::render(opt$reportTemplate,
                  output_file = name,
                  output_dir = file.path(opt$baseDir, opt$resultsDir, opt$prefix),
                  params=list(selfname = name,
                              ntop = opt$ntop,
                              pcutoff = opt$pcutoff,
                              lfcthreshold = opt$lfcthreshold,
                              baseDir = opt$baseDir,
                              resultsDir = file.path(opt$resultsDir, opt$prefix),
                              prefix = prefix,
                              derds = ddsrds,
                              minMean = opt$minMean,
                              minSingle = opt$minSingle,
                              filterBaseMean = opt$bmF,
                              reducedFormula = opt$reducedFormula,
                              longlabel=opt$label,
                              covars=covars,
                              widths = featLens,
                              specnorm=opt$specnorm,
                              rlog=!opt$vst,
                              comparisons=comparisons)
                  )





