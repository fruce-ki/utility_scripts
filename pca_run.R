#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec = matrix(c(
  'baseDir',        'b', 1, "character", "Full path to base directory of the project. All other paths relative from here.",
  'createID',       'C', 0, "logical",   "Non-summarized intron or exon counts by featureCounts. The first column is not unique IDs because genes are repeated. (FALSE)",
  'countsFile',     'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples.",
  'forVar',         'F', 1, "character", "Looping variable (ie. carry out the analysis separately for each level of this var).",
  'specnorm',       'G', 1, "character", "Special normalisation. Gene ids that match this pattern will not be included for calculation of library sizes for TPM/RPM normalisation.",
  'help',           'h', 0, "logical",   "Help",
  'idcol',          'i', 1, "integer",   "ID column to use (1). The others will be removed. See also -I and -F.",
  'nidcols',        'I', 1, "integer",   "Number of non-count columns at the start of the table (1). See also -i and -W.",
  'minMean',        'M', 1, "integer",   "Minimum mean count across all samples (10).",
  'minSingle',      'm', 1, "integer",   "Minimum count in at least one single sample (100).",
  'nhit',           'n', 1, "integer",   "Number of hits to report (10).",
  'resultsDir',     'o', 1, "character", "Directory in which to save the report (.).",
  'RDSoutdir'    ,  'r', 1, "character", "Directory in which to save the raw computed objects (.).",
  'samplesFile',    's', 1, "character", "Tab-separated table with `Sample` column followed by the variable columns.",
  'reportTemplate', 'T', 1, "character", "Full path to template Rmd file (~/utility_scripts/pca_report_template.Rmd).",
  'topVars',        'v', 1, "integer",   "How many genes to use, ranked by descending Coefficient of Variation. (500)",
  'widthsFile',     'w', 1, "character", "Full path to table of feature lengths. 2 columns: feature ID and length.  the IDs must match the count IDs and be unique.",
  'widthsCol',      'W', 1, "integer",   "Column in the counts file that contains feature lengths. This should be among the non-count columns, see also -I. It is an alternative to providing a separate lengths file."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(baseDir="/Volumes/groups/busslinger/Kimon/tanja/R14425_RNAseq", countsFile="process/featureCounts/intron_genecounts.txt", createID=FALSE, samplesFile="description/covars_pca.txt", resultsDir="results/PCA", idcol=1, nidcols=6, minMean=50, minSingle=100, nhit=15, reportTemplate="/Volumes/groups/busslinger/Kimon/tanja/R14425_RNAseq/code/pca_report_template.Rmd", widthsCol=6, specnorm='^ENSMUSG|^TCR|^IGH', topVars=500, forVar='Design')

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$reportTemplate) ) {
  opt$reportFile <- '~/utility_scripts/pca_report_template.Rmd'
}

stopifnot(!is.null(opt$baseDir))
stopifnot(!is.null(opt$countsFile))
stopifnot(!is.null(opt$samplesFile))

if (is.null(opt$resultsDir))  opt$resultsDir <- '.'
if (is.null(opt$createID))  opt$createID <- FALSE
if ((!is.null(opt$specnorm)) & opt$specnorm == "NULL")  opt$specnorm <- NULL
if (is.null(opt$nidcols)) opt$nidcols <- 1L
if (is.null(opt$idcol)) opt$idcol <- 1L
if (is.null(opt$minMean)) opt$minMean <- 10L
if (is.null(opt$minSingle)) opt$minSingle <- 100L
if (is.null(opt$nhit))  opt$nhit <- 10L
if (is.null(opt$topVars))  opt$topVars <- 500L
if ((!is.null(opt$forVar)) && opt$forVar == "NULL") opt$forVar <- NULL



dir.create(file.path(opt$baseDir, opt$resultsDir), recursive=TRUE)
if (!is.null(opt$RDSoutdir)) {
  opt$RDSdir <- file.path(opt$baseDir, opt$RDSoutdir)
  dir.create(opt$RDSoutdir, recursive=TRUE)
}


# Counts
counts <- fread(file.path(opt$baseDir, opt$countsFile), colClasses="character")  # prevent any numeric IDs converting a string column to integer
counts <- unique(counts)         # non-summarized featureCounts have repeated exons from sharing among transcripts.
xref <- counts[, 1:opt$nidcols]
samples <- names(counts)[c(opt$nidcols+1):length(counts)]
counts <- matrix(as.numeric(as.matrix(counts[, c(opt$nidcols+1):length(counts), with=FALSE])), ncol=length(counts) - opt$nidcols )
colnames(counts) <- samples  # because the type conversion lost them
counts[is.na(counts)] <- 0
if (opt$createID) {
  xref[, id := paste0(Geneid, '.', Start, '_', End)]
  if (any(duplicated(xref$id)))
    stop()
  rownames(counts) <- xref$id
} else {
  rownames(counts) <- xref[[opt$idcol]]
}

# Feature lengths
if (!is.null(opt$widthsCol)) {
  featLens <- data.table(ID=rownames(counts),
                         LEN=as.integer(xref[[opt$widthsCol]]))
  setkey(featLens, ID)
} else if (!is.null(opt$widthsFile)) {
  featLens <- fread(file.path(opt$baseDir, opt$widthsFile), sep="\t", header=TRUE, check.names = FALSE, colClasses=c("character", "numeric"))
  setkey(featLens, ID)
} else {
  featLens <- NULL
}


# TPM/RPM
normscale <- function(counts, featLens=NULL, specnorm=opt$specnorm){
  if (!is.null(featLens)) {
    # Reorder/subset
    setkeyv(featLens, names(featLens)[1])
    featLens <- featLens[rownames(counts)]
    stopifnot(all(!is.na(featLens[, 2])))  # incomplete lengths
    # Scale by feature size
    TPMs <- sweep(counts, 1, featLens[[2]], `/`)
  }
  # Scale by sequencing depth
  if (!is.null(specnorm)) {
    colsums <- colSums(TPMs[!grepl(specnorm, rownames(TPMs), ignore.case=TRUE, perl=TRUE), ], na.rm=TRUE)  # MOUSE only!
  } else {
    colsums <- colSums(TPMs, na.rm=TRUE)
  }
  TPMs <- sweep(TPMs, 2, colsums, `/`) * 1e6
  
  return(TPMs)
}

tpm <- normscale(counts, featLens)
TPM <- data.table(id=rownames(tpm), tpm)
widths <- !all(is.null(opt$widthsFile), is.null(opt$widthsCol))
setnames(TPM, c('id', paste0(colnames(tpm), ifelse(widths, '.TPM', '.RPM')) ))
fwrite(TPM, file=sub(".txt$|.tsv$", ".scaled.txt", file.path(opt$baseDir, opt$countsFile)), quote=FALSE, sep="\t")


# Covariates
covars <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", header=TRUE, check.names = FALSE, colClasses="character")
counts <- counts[, covars$Sample]    # Order and subset columns

# Repeat analysis for every level of the looping variable.
subcovars <- list()
if (!is.null(opt$forVar)) {
  # in case not all combinations of variables exist
  covars[[opt$forVar]] <- droplevels(as.factor(covars[[opt$forVar]]))  # sometimes it is a factor sometimes it isn,t
  
  subcovars <- lapply(levels(covars[[opt$forVar]]), function(V) {
    # V <- levels(covars[[opt$forVar]])[1]
    
    # Select relevant rows and drop the column
    coldata <- covars[covars[[opt$forVar]] == V, ]
    coldata[opt$forVar] <- NULL
    return(coldata)
  })
  names(subcovars) <- levels(covars[[opt$forVar]])
}
# Always also do the full data
subcovars <- c(list(covars), subcovars)
names(subcovars)[1] <- 'default'


for (V in names(subcovars)){
  # V <- names(subcovars)[1]
  
  # Fire up the Rmd report
  rmarkdown::render(opt$reportTemplate,
                    output_file = sub('.txt|.tsv', paste0('.', gsub('[^A-Za-z0-9]+', '_', V), '.pca.html'), basename(opt$countsFile)),
                    output_dir = file.path(opt$baseDir, opt$resultsDir),
                    params=list(cts = file.path(opt$baseDir, opt$countsFile),
                                tpm=tpm[, subcovars[[V]]$Sample],
                                # covars = file.path(opt$baseDir, opt$samplesFile),
                                covars = subcovars[[V]],
                                RDSdir = file.path(opt$baseDir, opt$RDSoutdir),
                                nidcols = opt$nidcols,
                                idcol = opt$idcol,
                                minMean = opt$minMean,
                                minSingle = opt$minSingle,
                                ntop = opt$nhit,
                                widths = widths,
                                createID=opt$createID,
                                specnorm=opt$specnorm,
                                topVars=opt$topVars,
                                loopVal=V)
  )
}
