#!/usr/bin/env Rscript

library(getopt)
library(DESeq2)

spec = matrix(c(
  'baseDir'      , 'b', 1, "character", "Base directory for everything (.)",
  'countsFile'   , 'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples",
  'help'         , 'h', 0, "logical",   "Help",
  'resultsDir'   , 'o', 1, "character", "Directory in which to save the contrast results. If omitted, only the RDS will be output.",
  'RDSoutdir'    , 'r', 1, "character", "Directory in which to save the raw DESeq2 object (./process)",
  'samplesFile'  , 's', 1, "character", "Tab-separated table with `sample` column followed by the variable columns",
  'control'      , 'x', 1, "character", "Comma separated value for each variable to use as reference for all comparisons, in the order variables are listed in samplesfile",
  'designFormula', 'd', 1, "character", "Design formula",
  'reducedFormula','R', 1, "character", "Reduced formula for LR test",
  'minCount',      'M', 1, "character", "Minimum number of reads across all samples combined for a gene to be considered (30)."
), byrow=TRUE, ncol=5)
opt = getopt(spec)
#opt <- list(baseDir='/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/M9179_quantseq', countsFile='process/counts/all_counts_named.tsv', resultsDir='./process/DE', RDSoutdir='./process/DE', samplesFile='./aux/CMcontrols.txt', control='post', designFormula='~ time', reducedFormula='~1')

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$baseDir    ) ) { 
  opt$baseDir    = '.'         
}
if ( is.null(opt$RDSoutdir ) ) { 
  opt$RDSoutdir = './process' 
}

if ( is.null(opt$minCount ) ) { 
  opt$minCount = 30L 
}

dir.create(file.path(opt$baseDir, opt$RDSoutdir))
if ( !is.null(opt$resultsDir ) ) { 
  dir.create(file.path(opt$baseDir, opt$resultsDir))
}

# Input
cts <- round(as.matrix(read.csv(file.path(opt$baseDir, opt$countsFile), sep="\t", header=TRUE, row.names=1)), digits=0)
coldata <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", header=TRUE, row.names='sample')

# Remove covariates with a single category (unnecessary since the design is provided explicitly).
# for (n in names(coldata)){
#   if ( length(unique(coldata[,n])) < 2 ) {
#     message(paste0("Dropped `", n, "` from the covariates table due to uniform value."))
#     coldata[[n]] <- NULL
#   }
# }

# Remove unused samples, and ensure the rest are in the correct order (simplifies subset comparisons from larger datasets).
cdn <- rownames(coldata)
ctsn <- colnames(cts)
if (!all(cdn %in% ctsn) ) {
  if( all(paste0('X',cdn) %in% ctsn) ) {  # Handle automatic X prepended to names starting with number
    cdn <- paste0('X',cdn)
  } else {
    print(cdn[!cdn %in% ctsn])
    stop("Error! Samples not in data!")
  }
}
cts <- cts[, match(cdn, ctsn)]

# Setup
vars = names(coldata)[ vapply(names(coldata), function(x) { grepl(x, opt$designFormula) }, logical(1)) ] # only variables used in the design
refs <- strsplit(opt$control, ',')[[1]]
names(refs) <- vars

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = as.formula(opt$designFormula))
# set reference condition(s)
for (n in vars){
  dds[[n]] <- relevel(dds[[n]], ref=refs[n])
}

autoname <- gsub(' |~', '', opt$designFormula)

# Prefilter out genes with very poor coverage.
# This is a speed and memory improvement. DESeq2 later internally applies stricter thresholds.
keep <- rowSums(counts(dds)) >= opt$minCount
dds <- dds[keep,]

# DE
if (!is.null(opt$reducedFormula)){
  dds <- DESeq(dds, test='LRT', reduced = as.formula(opt$reducedFormula))  # LR test
} else {
  dds <- DESeq(dds) # Wald test
}

saveRDS(dds, file=file.path(opt$baseDir, opt$RDSoutdir, paste0(autoname, '_deseq2.RDS')))

# Plain results
# res_e5 <- results(dds, contrast=c('condition', 'N1', 'E5'))
# res_e7 <- results(dds, contrast=c('condition', 'N1', 'E7'))
# res_m1 <- results(dds, contrast=c('condition', 'N1', 'M1'))

# Shrunk LFC supposedly better for plotting and ranking
if( !is.null(opt$resultsDir) ) {
  coefficients <- resultsNames(dds)
  for (name in coefficients[2:length(coefficients)]) {
    res <- lfcShrink(dds, coef=name, type='apeglm')       # apeglm is the type recommended by DESeq2
    
                                                          # Significance and counts are kept the same, but lfc is shrunk.
    res$mlog10p <- -log10(res$padj)
    write.table(data.frame(row_ID=rownames(res), res), file=file.path(opt$baseDir, opt$resultsDir, paste0(autoname, '_', name, '.tsv')),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
}
