#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'help',           'h', 0, "logical",   "Help",
  'countsFile',     'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples.",
  'samplesFile',    's', 1, "character", "Tab-separated table with `Sample` column followed by the variable columns.",
  'resultsDir',     'o', 1, "character", "Directory in which to save the report (.).",
  'RDSoutdir'    ,  'r', 1, "character", "Directory in which to save the raw computed objects (.).",
  'idcol',          'i', 1, "integer",   "ID column to use (1). The others will be removed. See also -I.",
  'nidcols',        'I', 1, "integer",   "Number of ID columns at the start of the table (1). See also -i.",
  'minMean',        'M', 1, "integer",   "Minimum mean count across all samples (10).",
  'minSingle',      'm', 1, "integer",   "Minimum count in at least one single sample (100).",
  'nhit',           'n', 1, "integer",   "Number of hits to report (10).",
  'reportTemplate', 'T', 1, "character", "Template Rmd file (~/utility_scripts/PCA_report_template.Rmd)."
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(reportTemplate='/groups/busslinger/Kimon/tanja/R12593_RNAseq/code/PCA_report_template.Rmd', countsFile='/groups/busslinger/Kimon/tanja/R12593_RNAseq/process/featureCounts/intron_se_genecounts.txt', samplesFile='/groups/busslinger/Kimon/tanja/R12593_RNAseq/description/covars_pca.txt', resultsDir='/groups/busslinger/Kimon/tanja/R12593_RNAseq/results/PCA', RDSoutdir='/groups/busslinger/Kimon/tanja/R12593_RNAseq/process/PCA', nidcols=6, idcol=1, nhit=25)


if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
 
if ( is.null(opt$reportTemplate) ) {
  opt$reportFile <- '~/utility_scripts/PCA_report_template.Rmd'
}

stopifnot(!is.null(opt$countsFile))
stopifnot(!is.null(opt$samplesFile))

if ( is.null(opt$resultsDir) ) {
  opt$resultsDir <- '.'
}
if ( is.null(opt$RDSoutdir ) ) { 
  opt$RDSoutdir <- '.' 
}

if ( is.null(opt$nidcols) ) {
  opt$nidcols <- 1L
}
if ( is.null(opt$idcol) ) {
  opt$idcol <- 1L
}

if ( is.null(opt$minMean) ) {
  opt$minMean <- 10L
}
if ( is.null(opt$minSingle) ) {
  opt$minSingle <- 100L
}

if ( is.null(opt$nhit) ){
  opt$nhit <- 10L
}


dir.create(opt$resultsDir, recursive=TRUE)
dir.create(opt$RDSoutdir, recursive=TRUE)

# Fire up the Rmd report
rmarkdown::render(opt$reportTemplate,
                  output_file = sub('.txt|.tsv', '_pca.html', basename(opt$countsFile)),
                  output_dir = opt$resultsDir,
                  params=list(cts = opt$countsFile,
                              covars = opt$samplesFile,
                              RDSdir = opt$RDSoutdir,
                              nidcols = opt$nidcols,
                              idcol = opt$idcol,
                              minMean = opt$minMean,
                              minSingle = opt$minSingle,
                              ntop = opt$nhit)
)
