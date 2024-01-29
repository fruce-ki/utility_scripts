#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec = matrix(c(
  'baseDir',        'b', 1, "character", "Full path to base directory of the project. All other paths relative from here.",
  'countsFile',     'f', 1, "character", "Salmon -like counts table with all the samples.",
  'excludeList',    'e', 1, "character", "Text file with vertical list of feature IDs to exclude from PCA.",
  'help',           'h', 0, "logical",   "Help",
  'idcol',          'i', 1, "integer",   "ID column to use (1). The others will be removed. See also -I and -F.",
  'nidcols',        'I', 1, "integer",   "Number of non-count columns at the start of the table (1). See also -i and -W.",
  'minMean',        'M', 1, "integer",   "Minimum mean count across all samples (10).",
  # 'minSingle',      'm', 1, "integer",   "Minimum count in at least one single sample (100).",
  'resultsDir',     'o', 1, "character", "Directory in which to save the report, relative to baseDir (.).",
  'samplesFile',    's', 1, "character", "Tab-separated table with `sample` column followed by the variable columns.",
  'sortVar',        'S', 1, "character", " a variable by which to order the samples. Otherwise the samplesFile order will be used.",
  'forVar',         'F', 1, "character", "Looping variable (ie. carry out the analysis separately for each level of this var).",
  'reportTemplate', 'T', 1, "character", "Full path to template Rmd file (~/utility_scripts/pca_report_template.Rmd).",
  'nhit',           'n', 1, "integer",   "Number of hits to report (10).",
  'topVars',        'v', 1, "integer",   "How many genes to use, ranked by descending Coefficient of Variation. (500)",
  'allSamples',     'a', 0, "logical",   "Do not include all the samples, include only the ones marked by 'use' (FALSE)."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(baseDir="/SCRATCH/PP2023011_SLC13A5", countsFile="tesslinz.salmon.merged.gene_tpm.tsv", samplesFile="samplesheet_tesslinz.differentialabundance.csv", resultsDir="tesslinz_diffex/PCA", allSamples = TRUE, sortVar = 'time', idcol=1, nidcols=2, minMean=5)
# opt <- list(baseDir="C:/Users/jack_/Downloads/", countsFile="tesslinz.salmon.merged.gene_tpm.tsv", samplesFile="samplesheet_tesslinz.differentialabundance.csv", resultsDir="tesslinz_diffex/PCA", sortVar = 'time', idcol=1, nidcols=2, minMean=5, reportTemplate="D:/Documents/GitHub/utility_scripts/pca_report_template.Rmd")

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$reportTemplate) ) {
  opt$reportTemplate <- '~/utility_scripts/pca_report_template.Rmd'
}

stopifnot(!is.null(opt$baseDir))
stopifnot(!is.null(opt$countsFile))
stopifnot(!is.null(opt$samplesFile))

if (is.null(opt$resultsDir))  opt$resultsDir <- '.'
if (is.null(opt$nidcols)) opt$nidcols <- 1L
if (is.null(opt$idcol)) opt$idcol <- 1L
if (is.null(opt$minMean)) opt$minMean <- 10L
if (is.null(opt$minSingle)) opt$minSingle <- 100L
if (is.null(opt$nhit))  opt$nhit <- 10L
if (is.null(opt$topVars))  opt$topVars <- 500L
if ((!is.null(opt$forVar)) && opt$forVar == "NULL") opt$forVar <- NULL
if ((!is.null(opt$groupVar)) && opt$groupVar == "NULL") opt$groupVar <- NULL
if(is.null(opt$allSamples)) opt$allSamples <- FALSE


dir.create(file.path(opt$baseDir, opt$resultsDir), recursive=TRUE)


# Genes not to include in PCA. For example, cell cycle.
exclusion <- NULL
if (!is.null(opt$excludeList))
  exclusion <- read.csv(opt$excludeList, header=FALSE)[[1]]

# Covariates
covars <- fread(file.path(opt$baseDir, opt$samplesFile), sep=",", header=TRUE, check.names = FALSE, colClasses="character")


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
  print(V)

  # Fire up the Rmd report
  rmarkdown::render(opt$reportTemplate,
                    output_file = sub('.txt|.tsv', paste0('.', gsub('[^A-Za-z0-9]+', '_', V), '.pca.html'), opt$countsFile),
                    output_dir = file.path(opt$baseDir, opt$resultsDir),
                    params = list(name = file.path(opt$baseDir, opt$countsFile),
                                  tpm = file.path(opt$baseDir, opt$countsFile),
                                  covars = subcovars[[V]],
                                  outdir = file.path(opt$baseDir, opt$resultsDir),
                                  nidcols = opt$nidcols,
                                  idcol = opt$idcol,
                                  minMean = opt$minMean,
                                  # minSingle = opt$minSingle,
                                  ntop = opt$nhit,
                                  topVars=opt$topVars,
                                  loopVal=V,
                                  excluded=exclusion,
                                  groupvar=opt$sortVar,
                                  onlyUsable=opt$allSamples)
  )
}

