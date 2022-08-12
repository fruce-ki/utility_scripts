#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'baseDir',        'b', 1, "character", "Full path to base directory of the project. All other paths relative from here.",
  'createID',       'C', 0, "logical",   "Non-summarized intron or exon counts by featureCounts. The first column is not unique IDs because genes are repeated. (FALSE)",
  'countsFile',     'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples.",
  'forvar',         'F', 1, "character", "Looping variable (ie. carry out the analysis separately for each level of this var).",
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

# opt <- list(baseDir="/Volumes/groups/busslinger/Kimon/anna/R13870_RNAseq", countsFile="process/featureCounts/introns_genecounts.txt", createID=FALSE, samplesFile="description/covars.txt", resultsDir="results/PCA", RDSoutdir="process/PCA", idcol=1, nidcols=6, minMean=50, minSingle=100, nhit=15, reportTemplate="/Volumes/groups/busslinger/Kimon/anna/R13870_RNAseq/code/pca_report_template.Rmd", widthsCol=6, specnorm='^ENSMUSG|^TCR|^IGH', topVars=500, forvar='NULL')

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
if (is.null(opt$RDSoutdir)) opt$RDSoutdir <- '.'
if (is.null(opt$createID))  opt$createID <- FALSE
if ((!is.null(opt$specnorm)) & opt$specnorm == "NULL")  opt$specnorm <- NULL
if (is.null(opt$nidcols)) opt$nidcols <- 1L
if (is.null(opt$idcol)) opt$idcol <- 1L
if (is.null(opt$minMean)) opt$minMean <- 10L
if (is.null(opt$minSingle)) opt$minSingle <- 100L
if (is.null(opt$nhit))  opt$nhit <- 10L
if (is.null(opt$topVars))  opt$topVars <- 500L
if ((!is.null(opt$forVar)) && opt$forvar == "NULL") opt$forVar <- NULL



dir.create(file.path(opt$baseDir, opt$resultsDir), recursive=TRUE)
dir.create(file.path(opt$baseDir, opt$RDSoutdir), recursive=TRUE)


# Covariates
covars <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", header=TRUE, check.names = FALSE, colClasses="character")

# Repeat analysis for every level of the looping variable.
subcovars <- list()
if (!is.null(opt$forar)) {
  # in case not all combinations of variables exist
  covars[[opt$forvar]] <- droplevels(as.factor(covars[[opt$forvar]]))  # sometimes it is a factor sometimes it isn,t
  
  subcovars <- lapply(levels(covars[[opt$forvar]]), function(V) {
    # V <- levels(covars[[opt$forvar]])[1]
    
    # Select relevant rows and drop the column
    coldata <- covars[covars[[opt$forvar]] == V, ]
    coldata[opt$forvar] <- NULL
    return(coldata)
  })
  names(subcovars) <- levels(covars[[opt$forvar]])
}
# Always also do the full data
subcovars <- c(subcovars, list(covars))
names(subcovars)[length(subcovars)] <- 'default'

for (V in names(subcovars)){
  # V <- names(subcovars)[1]
  coldata <- subcovars[[V]]
  
  # Fire up the Rmd report
  rmarkdown::render(opt$reportTemplate,
                    output_file = sub('.txt|.tsv', paste0('.', gsub('[^A-Za-z0-9]+', '_', V), '.pca.html'), basename(opt$countsFile)),
                    output_dir = file.path(opt$baseDir, opt$resultsDir),
                    params=list(cts = file.path(opt$baseDir, opt$countsFile),
                                # covars = file.path(opt$baseDir, opt$samplesFile),
                                covars = coldata,
                                RDSdir = file.path(opt$baseDir, opt$RDSoutdir),
                                nidcols = opt$nidcols,
                                idcol = opt$idcol,
                                minMean = opt$minMean,
                                minSingle = opt$minSingle,
                                ntop = opt$nhit,
                                widthsFile = opt$widthsFile,
                                widthsCol=opt$widthsCol,
                                createID=opt$createID,
                                specnorm=opt$specnorm,
                                topVars=opt$topVars,
                                loopVal=V)
  )
}
