#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'help'         , 'h', 0, "logical",   "Help",
  'countsFile'   , 'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples",
  'samplesFile'  , 's', 1, "character", "Tab-separated table with `sample` column followed by the variable columns",
  'nidcols',       'I', 1, "numeric",   "Number of ID columns at the start of the table (1). see also -i",
  'idcol',         'i', 1, "numeric",   "ID column to use (1). The others will be removed. See also -I.",
  'resultsDir'   , 'o', 1, "character", "Directory in which to save the contrast results. If omitted, only the RDS will be output.",
  'reportTemplate','T', 1, "character", "Template Rmd file for DE report."
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(countsFile='/Volumes/groups/zuber/zubarchive/USERS/Kimon/anja/M9186_quantseq/process_quant/all_counts_rpm_xref.txt', samplesFile='/Volumes/groups/zuber/zubarchive/USERS/Kimon/anja/M9186_quantseq/description/covars.txt', resultsDir='/Volumes/groups/zuber/zubarchive/USERS/Kimon/anja/M9186_quantseq/results_quant/PCA', nidcols=3, idcol=2)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$reportFile) ) {
  opt$reportFile <- '~/utility_scripts/PCA_report_template.Rmd'
}

if ( is.null(opt$nidcols ) ) { 
  opt$nidcols <- 1
}

if ( is.null(opt$resultsDir    ) ) { 
  opt$resultsDir <- '.'         
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
                              pdf = file.path(opt$resultsDir, sub('.txt|.tsv', '_pca.pdf', basename(opt$countsFile))),
                              prefix = file.path(opt$resultsDir, sub('.txt|.tsv', '', basename(opt$countsFile))) )
)
