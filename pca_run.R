#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'help',           'h', 0, "logical",   "Help",
  'countsFile',     'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples.",
  'samplesFile',    's', 1, "character", "Tab-separated table with `sample` column followed by the variable columns.",
  'nvar',           'v', 1, 'integer',   "Number of most variable features to use (500).",
  'nidcols',        'I', 1, "integer",   "Number of ID columns at the start of the table (1). See also -i.",
  'idcol',          'i', 1, "integer",   "ID column to use (1). The others will be removed. See also -I.",
  'resultsDir',     'o', 1, "character", "Directory in which to save the contrast results. If omitted, only the RDS will be output.",
  'reportTemplate', 'T', 1, "character", "Template Rmd file for DE report."
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(countsFile='/groups/zuber/zubarchive/USERS/Kimon/markus/OTI_vivo_pdac/process/all_counts_rpm_xref.txt', samplesFile='/groups/zuber/zubarchive/USERS/Kimon/markus/OTI_vivo_pdac/description/covars_only-no.txt', resultsDir='/groups/zuber/zubarchive/USERS/Kimon/markus/OTI_vivo_pdac/results/PCA_only-no_test', nidcols=3, idcol=2, nvar=500)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$reportFile) ) {
  opt$reportFile <- '~/utility_scripts/PCA_report_template.Rmd'
}

if ( is.null(opt$nidcols) ) {
  opt$nidcols <- 1L
}
if ( is.null(opt$idcol) ) {
  opt$idcol <- 1L
}

if ( is.null(opt$resultsDir) ) {
  opt$resultsDir <- '.'
}

if ( is.null(opt$nvar) ) {
  opt$nvar <- 500L
}


dir.create(opt$resultsDir, recursive=TRUE)

# Fire up the Rmd report
rmarkdown::render(opt$reportFile,
                  output_file = sub('.txt|.tsv', '_pca.html', basename(opt$countsFile)),
                  output_dir = opt$resultsDir,
                  params=list(cts = opt$countsFile,
                              covars = opt$samplesFile,
                              nidcols = opt$nidcols,
                              idcol = opt$idcol,
                              topvars = opt$nvar,
                              pdf = file.path(opt$resultsDir, sub('.txt|.tsv', '_pca.pdf', basename(opt$countsFile))),
                              prefix = file.path(opt$resultsDir, sub('.txt|.tsv', '', basename(opt$countsFile))) )
)
