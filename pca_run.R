#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec = matrix(c(
  'baseDir',        'b', 1, "character", "Full path to base directory of the project. All other paths relative from here.",
  'createID',       'C', 0, "logical",   "Non-summarized intron or exon counts by featureCounts. The first column is not unique IDs because genes are repeated. (FALSE)",
  'countsFile',     'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples.",
  'forVar',         'F', 1, "character", "Looping variable (ie. carry out the analysis separately for each level of this var).",
  'excludeList',    'g', 1, "character", "Text file with vertical list of feature IDs to exclude from PCA (but not from TPM calculation).",
  'specnorm',       'G', 1, "character", "Special normalisation. Feature IDs that match this pattern will not be included for calculation of library sizes for TPM/RPM normalisation, and will also not be included in the PCA.",
  'help',           'h', 0, "logical",   "Help",
  'idcol',          'i', 1, "integer",   "ID column to use (1). The others will be removed. See also -I and -F.",
  'nidcols',        'I', 1, "integer",   "Number of non-count columns at the start of the table (1). See also -i and -W.",
  'minMean',        'M', 1, "integer",   "Minimum mean count across all samples (10).",
  'minSingle',      'm', 1, "integer",   "Minimum count in at least one single sample (100).",
  'nhit',           'n', 1, "integer",   "Number of hits to report (10).",
  'resultsDir',     'o', 1, "character", "Directory in which to save the report (.).",
  'samplesFile',    's', 1, "character", "Tab-separated table with `sample` column followed by the variable columns.",
  'reportTemplate', 'T', 1, "character", "Full path to template Rmd file (~/utility_scripts/pca_report_template.Rmd).",
  'topVars',        'v', 1, "integer",   "How many genes to use, ranked by descending Coefficient of Variation. (500)",
  'widthsFile',     'w', 1, "character", "Full path to table of feature lengths. 2 columns: feature ID and length.  the IDs must match the count IDs and be unique.",
  'widthsCol',      'W', 1, "integer",   "Column in the counts file that contains feature lengths. This should be among the non-count columns, see also -I. It is an alternative to providing a separate lengths file."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(baseDir="/SCRATCH/PP2023011_SLC13A5", countsFile="tesslinz.salmon.merged.gene_tpm.tsv", createID=FALSE, samplesFile="samplesheet_tesslinz.differentialabundance.csv", resultsDir=".", RDSoutdir='.', idcol=1, nidcols=2, minMean=5, minSingle=20, reportTemplate="~/utility_scripts/pca_report_template.Rmd")

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
if ((!is.null(opt$specnorm)) && opt$specnorm == "NULL")  opt$specnorm <- NULL
if (is.null(opt$nidcols)) opt$nidcols <- 1L
if (is.null(opt$idcol)) opt$idcol <- 1L
if (is.null(opt$minMean)) opt$minMean <- 10L
if (is.null(opt$minSingle)) opt$minSingle <- 100L
if (is.null(opt$nhit))  opt$nhit <- 10L
if (is.null(opt$topVars))  opt$topVars <- 500L
if ((!is.null(opt$forVar)) && opt$forVar == "NULL") opt$forVar <- NULL


dir.create(file.path(opt$baseDir, opt$resultsDir), recursive=TRUE)


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
normscale <- function(counts, featLens = NULL, specnorm = NULL){
  # counts <- M
  TPMs <- NULL
  if (!is.null(featLens)) {
    stopifnot(all(is.finite(featLens)))
    stopifnot(length(featLens) == nrow(counts))
    # Scale by feature size
    message("Scale to TPM.")
    TPMs <- sweep(counts, 1, featLens, `/`)
  } else {
    message("Scale to RPM.")
    TPMs <- counts
  }
  # Scale by sequencing depth
  if (!is.null(specnorm)) {
    message("Scale with special exclusions.")
    colsums <- colSums(TPMs[!grepl(specnorm, rownames(TPMs), ignore.case=TRUE, perl=TRUE), ], na.rm=TRUE)
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
#fwrite(TPM, file=sub(".txt$|.tsv$", ".scaled.txt", file.path(opt$baseDir, opt$countsFile)), quote=FALSE, sep="\t")


# Genes not to include in PCA. For example, cell cycle.
exclusion <- NULL
if (!is.null(opt$excludeList))
  exclusion <- read.csv(opt$excludeList, header=FALSE)[[1]]


# Covariates
covars <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep=",", header=TRUE, check.names = FALSE, colClasses="character")
counts <- counts[, covars$sample]    # Order and subset columns

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
                                tpm=tpm[, subcovars[[V]]$sample],
                                counts = counts,
                                covars = subcovars[[V]],
                                outdir = file.path(opt$baseDir, opt$resultsDir),
                                nidcols = opt$nidcols,
                                idcol = opt$idcol,
                                minMean = opt$minMean,
                                minSingle = opt$minSingle,
                                ntop = opt$nhit,
                                widths = widths,
                                createID=opt$createID,
                                specnorm=opt$specnorm,
                                topVars=opt$topVars,
                                loopVal=V,
                                excluded=exclusion)
  )
}
