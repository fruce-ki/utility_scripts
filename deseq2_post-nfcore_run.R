#!/usr/bin/env Rscript

###
# Not intended to be given to nf-core/differentialabundance directly to replace the default template.
# Instead, meant to be run aseparately after completion of the pipeline. It ralies on the outputs of the pipeline
###

library(getopt)

spec = matrix(c(
  'baseDir',        'b', 1, "character", "Full path to base directory of the nf-core/differentialabundance results. All other paths relative from here.",
  'countsFile',     'f', 1, "character", "Salmon-formatted gene counts file.",
  'tpmFile',        'F', 1, "character", "Salmon-formatted gene TPMs file.",
  'samplesFile',    's', 1, "character", "Tab-separated table with `sample` column followed by the variable columns.",
  'contrasts',      'c', 1, "character", "File with the contrasts for nf-core/differentialabundance.",
  'genelist',       'g', 1, "character", "File with (vertical) list of genes to highlight, instead of top DE genes.",
  'minCount',       'm', 1, "integer",   "Minimum required mean count in at least one condition (50).",
  'minTPM',         'M', 1, "numeric",   "Minimum required mean TPM in at least one single sample (5).",
  'minLFC',         'l', 1, "numeric",   "Minimum required log2 fold-change (0.5).",
  'maxQ',           'p', 1, "numeric",   "FDR level (0.05).",
  'nhit',           'n', 1, "integer",   "Number of hits to report (50).",
  'pad',            'z', 1, "logical",   "Disable padding of zeros.",
  'suffix',         'x', 1, "character", "A custom string to add to the report file name.",
  'reportTemplate', 'T', 1, "character", "Full path to template Rmd file (~/utility_scripts/deseq2_post-nfcore_report_template.Rmd).",
  'help',           'h', 0, "logical",   "Help."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


if (is.null(opt$minCount))
  opt$minCount <- 50
if (is.null(opt$minTPM))
  opt$minTPM <- 5
if (is.null(opt$minLFC))
  opt$minLFC <- 0.5
if (is.null(opt$maxQ))
  opt$maxQ <- 0.05
if (is.null(opt$nhit))
  opt$nhit <- 50
if (is.null(opt$reportTemplate))
  opt$reportTemplate <- "~/utility_scripts/deseq2_post-nfcore_report_template.Rmd"
if (is.null(opt$pad))
  opt$pad <- FALSE


contrasts <- read.table(opt$contrasts, sep = ',', header = TRUE)

for (comparison in contrasts$id) {
  rmarkdown::render(opt$reportTemplate,
                    output_file = paste0(comparison, opt$suffix, '.deseq2_fgsea.report.html'),
                    output_dir = file.path(opt$baseDir, 'report'),
                    params=list(comparison = comparison,
                                da_path = opt$baseDir,
                                counts = opt$countsFile,
                                tpm = opt$tpmFile,
                                covars = opt$samplesFile,
                                highlight_n = opt$nhit,
                                min_count = opt$minCount,
                                min_tpm = opt$minTPM,
                                min_lfc = opt$minLFC,
                                max_q = opt$maxQ,
                                genelist = opt$genelist,
                                pad = !opt$pad
                    )  )
}