#!/users/kimon.froussios/miniconda3/envs/bioinfo/bin/Rscript

library(getopt)

spec = matrix(c(
  'baseDir'      , 'b', 2, "character", "Base directory for everything (.)",
  'countsFile'   , 'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples",
  'help'         , 'h', 0, "logical",   "Help",
  'resultsDir'   , 'o', 2, "character", "Directory in which to save the contrast results (./results)",
  'RDSoutdir'    , 'r', 2, "character", "Directory in which to save the raw DESeq2 object (./process)",
  'samplesFile'  , 's', 1, "character", "Tab-separated table with `sample` and `condition` (2 columns)",
  'control'      , 'x', 1, "character", "Value of condition to use as reference for all comparisons"
), byrow=TRUE, ncol=5)
opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$baseDir    ) ) { opt$baseDir    = '.'         }
if ( is.null(opt$resultsDir ) ) { opt$resultsDir = './results' }
if ( is.null(opt$RDSoutdirir ) ) { opt$RDSoutdir = './process' }

# Input
cts <- round(as.matrix(read.csv(file.path(opt$baseDir, opt$countsFile), sep="\t", header=TRUE, row.names='row_ID')), digits=0)
coldata <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", header=TRUE, row.names='sample')

library(DESeq2)

# Setup
dds <- DESeqDataSetFromMatrix(countData = cts,
                               colData = coldata,
                               design = ~ condition)
dds$condition <- relevel(dds$condition, ref=opt$control)

autoname <- paste(levels(dds$condition), collapse='-')

# Prefilter
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]

# DE
dds <- DESeq(dds)
saveRDS(dds, file=file.path(opt$baseDir, opt$RDSoutdir, paste0(autoname, '_deseq2.RDS')))

# Plain results
# res_e5 <- results(dds, contrast=c('condition', 'N1', 'E5'))
# res_e7 <- results(dds, contrast=c('condition', 'N1', 'E7'))
# res_m1 <- results(dds, contrast=c('condition', 'N1', 'M1'))

# Shrunk LFC supposedly better for plotting and ranking
coefficients <- resultsNames(dds)
for (name in coefficients[2:length(coefficients)]) {
  res <- lfcShrink(dds, coef=name, type='apeglm')
  write.csv(res, file=file.path(opt$baseDir, opt$resultsDir, paste0(autoname, '_', name, '.tsv')),
            sep="\t", row.names=TRUE, col.names=TRUE)
}
