#!/usr/bin/env Rscript

library(getopt)
library(DESeq2)
library(data.table)

spec = matrix(c(
  'all',           'A', 0, "logical",   "Enable round-robin-style pairwise comparisons. Otherwise comparisons only against reference. (False)", 
  'bmF',           'B', 0, "logical",   "Disable baseMean filtering. (False)",
  'baseDir',       'b', 1, "character", "Base directory for everything (.)",
  'lfcthreshold',  'c', 1, "numeric",   "Log2 fold-change threshold (1)",
  'createID',      'C', 0, "logical",   "Non-summarized intron or exon counts by featureCounts. The first column is not unique IDs because genes are repeated.",
  'designFormula', 'd', 1, "character", "Design formula",
  'countsFile',    'f', 1, "character", "Tab-separated table of raw counts with `row_ID` and all the samples.",
  'forvar',        'F', 1, "character", "Looping variable (ie. carry out the design formula for each level of this var, separately).",
  'specnorm',       'G', 1, "character", "Special normalisation. Gene ids that match this pattern will not be included for calculation of library sizes for TPM/RPM normalisation.",
  'idcol',         'i', 1, "numeric",   "ID column to use (1). The others will be removed. See also -I.",
  'nidcols',       'I', 1, "numeric",   "Number of ID columns at the start of the table (1). See also -i.",
  'prescaled',     'k', 0, "logical",   "Don't let Deseq2 rescale the libraries. (False)",
  'label',         'l', 0, "logical",   "Add comparison details to the standard column names. Useful for merging multiple outputs.",
  'sellev',        'L', 1, "character", "Selection variable level. Used together with -S, to reduce the samplesFile to just the rows having that value in that variable.",
  'ntop',          'n', 1, "numeric",   "Number of hits to highlight (50)",
  'minMean',       'M', 1, "numeric",   "Minimum mean number of reads across all samples combined for a gene to be considered (10).",
  'minSingle',     'm', 1, "numeric",   "Minimum number of reads in any sinlge sample for a gene to be considered if the minMean is not met (100).",
  'resultsDir',    'o', 1, "character", "Directory in which to save the contrast results. If omitted, only the RDS will be output.",
  'pcutoff',       'p', 1, "numeric",   "P-value cutoff (0.05)",
  'prefix',        'P', 1, "character", "Prefix to add to output file name and column names (''). ",
  'RDSoutdir',     'r', 1, "character", "Directory in which to save the raw DESeq2 object (./process)",
  'reducedFormula','R', 1, "character", "Reduced formula for LR test",
  'samplesFile',   's', 1, "character", "Tab-separated table with `Sample` column followed by the variable columns. Values must NOT contain punctuation other than '_' !!!",
  'selvar',        'S', 1, "character", "Selection variable (always together with sellev). It is used together with -L. It is applied before -F.",
  'reportTemplate','T', 1, "character", "Template Rmd file for DE report.",
  'vst',           'v', 0, "logical",   "Plot VST counts instead of rlog2",
  'widthsFile',    'w', 1, "character", "Full path to table of feature lengths",
  'control',       'x', 1, "character", "Comma separated value for each variable to use as reference for all comparisons, in the order variables are listed in samplesfile",
  'widthsCol',      'W', 1, "integer",   "Column in the counts file that contains feature lengths. This should be among the non-count columns, see also -I. It is an alternative to providing a separate lengths file."
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(createID=FALSE, baseDir='/Volumes/groups/busslinger/Kimon/tanja/R12593_RNAseq', countsFile='process/featureCounts/intron_genecounts.txt', resultsDir='results/DE', RDSoutdir='process/DE', samplesFile='description/covars.txt', control='noTir1,0h,20210924,Tir1_end,2_5h', designFormula='~Condition', forvar="Group", selvar=NULL, sellev=NULL, minMean=50, minSingle=100, lfcthreshold=1, nidcols=6, idcol=1, ntop=50, bmF=FALSE, pcutoff=0.05, all=TRUE, prescaled=FALSE, prefix='intron_genecounts', label=TRUE, widthsCol=6, reportTemplate="/Volumes/groups/busslinger/Kimon/tanja/R12593_RNAseq/code/deseq2_report_template.Rmd")

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

stopifnot(!is.null(opt$countsFile))
stopifnot(!is.null(opt$samplesFile))

if (is.null(opt$reportTemplate))  opt$reportTemplate <- '~/utility_scripts/deseq2_report_template.Rmd'
if (is.null(opt$createID))  opt$createID <- FALSE 
if ((!is.null(opt$specnorm)) && opt$specnorm == "NULL") opt$specnorm <- NULL
if (is.null(opt$prescaled)) opt$prescaled <- FALSE
if (is.null(opt$label)) opt$label <- FALSE
if (is.null(opt$baseDir)) opt$baseDir <- '.'
if (is.null(opt$RDSoutdir)) opt$RDSoutdir <- './process' 
if (is.null(opt$resultsDir))  opt$resultsDir <- './results' 
if (is.null(opt$minMean)) opt$minMean <- 10
if (is.null(opt$minSingle)) opt$minSingle <- 100
if (is.null(opt$ntop))  opt$ntop <- 50
if (is.null(opt$pcutoff)) opt$pcutoff <- 0.05
if (is.null(opt$lfcthreshold))  opt$lfcthreshold <- 1
if (is.null(opt$nidcols)) opt$nidcols <- 1
if (is.null(opt$idcol)) opt$idcol <- 1
if (is.null(opt$bmF)) opt$bmF <- FALSE
opt$bmF <- !opt$bmF    # Convert -B flag functionality from switch-on to switch-off. 
if (is.null(opt$all)) opt$all <- FALSE
if (is.null(opt$vst)) opt$vst <- FALSE

# Output destinations
dir.create(file.path(opt$baseDir, opt$RDSoutdir, opt$prefix), recursive=TRUE)
dir.create(file.path(opt$baseDir, opt$resultsDir, opt$prefix), recursive=TRUE)


# Covariates
covars <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", row.names="Sample", header=TRUE, check.names = FALSE, colClasses="character")

# Counts
cts <- fread(file.path(opt$baseDir, opt$countsFile), sep="\t", header=TRUE, check.names = FALSE, colClasses="character")
# Setting all columns to character prevents numeric IDs at the top of the ID columns from auto-defining these columns as numeric when there could be character IDs further down that would otherwise become NaN.
cts <- unique(cts)         # non-sumxmarized featureCounts have repeated exons from sharing among transcripts.
xref <- cts[, 1:opt$nidcols]
cts <- as.matrix(cts[, rownames(covars), with=FALSE])    # Order and subset columns
cts <- matrix(as.numeric(cts), ncol=ncol(cts))
colnames(cts) <- row.names(covars)  # because the type conversion lost them
cts[is.na(cts)] <- 0
if (opt$createID) {
  xref[, id := paste0(Geneid, '.', Start, '_', End)]
  if (any(duplicated(xref$id)))
    stop()
  rownames(cts) <- xref$id
} else {
  rownames(cts) <- xref[[opt$idcol]]
}

# Feature lengths
if (!is.null(opt$widthsCol)) {
  featLens <- data.table(ID=rownames(cts),
                         LEN=as.integer(xref[[opt$widthsCol]]))
} else if (!is.null(opt$widthsFile)) {
  featLens <- fread(opt$widthsFile, sep="\t", header=TRUE, check.names = FALSE, colClasses=c("character", "numeric"))
  # Match order of feature lengths to order of features
  setkeyv(featLens, names(featLens)[1])
  featLens <- featLens[rownames(cts), ]
} else {
  featLens <- NULL
}

stopifnot(all(rownames(cts) == featLens[[1]]))
      

## DE prep stuff begins here. Done with the imports.

    
# Identify control level for the variables (it's less command editing to always provide all the controls, even the ones for variables that will be dropped).
refs <- strsplit(opt$control, ',')[[1]]
names(refs) <- names(covars)

# Apply sample selection
if ( !is.null(opt$selvar) ) {
  stopifnot(!is.null(opt$sellev) & opt$sellev %in% covars[[opt$selvar]])
  
  # Select the samples and then drop the variable, since it has only one value left after selection.
  covars <- covars[covars[[opt$selvar]] == opt$sellev, ]
  covars[opt$selvar] <- NULL
}


#################
# Repeat analysis for every level of the looping variable.
subcovars <- NULL
if (!is.null(opt$forvar)) {
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
} else {
  subcovars <- list(covars)
  names(subcovars) <- 'default'
}


for (V in names(subcovars)){
  # V <- names(subcovars)[1]
  coldata <- subcovars[[V]]
  
  # Identify the (remaining) variables.
  vars <- names(coldata)
  # Identify the ones relevant to the formula(s), and drop the rest.
  relevant <- vapply(names(coldata), function(x) { grepl(x, opt$designFormula) }, logical(1))
  vars <- names(coldata)[relevant]
  refs <- refs[vars]
  coldata <- coldata[, c(relevant), drop=FALSE]
  
  # Remove unused samples from the counts table, 
  # and ensure the rest are in the same order as the covariates table.
  cdn <- rownames(coldata)
  ctsn <- colnames(cts)
  subcts <- cts[, ctsn[match(cdn, ctsn)]]
  subcts[is.na(subcts)] <- 0
  
  # DESeq2 begins here
  DDS <- DESeqDataSetFromMatrix(countData = round(subcts, digits=0),
                                colData = coldata,
                                design = as.formula(opt$designFormula))
  # Set reference condition(s)
  for (n in vars){
    DDS[[n]] <- relevel(DDS[[n]], ref=refs[n])
  }
  if (opt$prescaled) {
    factors <- rep(1, times=length(cdn))
    names(factors) <- cdn
    sizeFactors(DDS) <- factors
  }
  
  # Prefilter out genes with very poor coverage.
  # DESeq2 internally applies additional thresholds.
  keep <- rowMeans(counts(DDS)) >= opt$minMean | rowSums(counts(DDS) >= opt$minSingle) >= 1
  DDS <- DDS[keep, ]
  # subrpm <- subrpm[keep, ]
  # subxref <- xref[keep, ]
  
  # Prepare file name.
  autoname <- NULL
  if (!is.null(opt$selvar)) {
    autoname <- paste(opt$selvar, opt$sellev, sep="_")
  }
  if (!is.null(opt$forvar)) {
    if (! is.null(autoname)){
      autoname <- paste(autoname, paste(opt$forvar, V, sep="_"), sep=".")
    } else {
      autoname <- paste(opt$forvar, V, sep="_")
    }
  }
  
  if (is.null(autoname) && is.null(opt$prefix)) {
    prefix <- gsub('[ ~]', '', opt$designFormula)
  } else if (!is.null(autoname) && is.null(opt$prefix)) {
    prefix <- paste(autoname, gsub(' ', '', opt$designFormula), sep='.')
  } else if (is.null(autoname) && !is.null(opt$prefix)) { 
    prefix <- paste(opt$prefix, gsub(' ', '', opt$designFormula), sep='.')
  } else {
    prefix <- paste(opt$prefix, gsub(' ', '', opt$designFormula), autoname, sep='.')
  }
  
  ddsrds <- file.path(opt$baseDir, opt$RDSoutdir, opt$prefix, paste0(prefix, '.deseq2data.RDS'))

  saveRDS(DDS, file=ddsrds)

  # Fire up an Rmd report
  name <- paste0(prefix, '.deseq2.html')
  rmarkdown::render(opt$reportTemplate, 
                    output_file = name,
                    output_dir = file.path(opt$baseDir, opt$resultsDir, opt$prefix),
                    params=list(selfname = name,
                                ntop = opt$ntop,
                                pcutoff = opt$pcutoff,
                                lfcthreshold = opt$lfcthreshold,
                                baseDir = opt$baseDir,
                                resultsDir = file.path(opt$resultsDir, opt$prefix),
                                prefix = prefix,
                                derds = ddsrds,
                                minMean = opt$minMean,
                                minSingle = opt$minSingle,
                                filterBaseMean = opt$bmF,
                                reducedFormula = opt$reducedFormula,
                                roundrobin=opt$all,
                                longlabel=opt$label,
                                covars=subcovars[[V]],
                                widths = featLens,
                                specnorm=opt$specnorm,
                                rlog=!opt$vst)
                    )
  
  # # DE
  # if (!is.null(opt$reducedFormula)){
  #   DDS <- DESeq(DDS, test='LRT', reduced = as.formula(opt$reducedFormula))  # LR test
  # } else {
  #   DDS <- DESeq(DDS) # Wald test
  # }
  #
  # # Compute contrasts with just the basic table outputs, no figures.
  # coefficients <- resultsNames(DDS)
  # if( !is.null(opt$resultsDir) ) {
  #   # Plain results
  #   for (name in coefficients[2:length(coefficients)]) {
  #     # name <- coefficients[2]
  #     
  #     res <- results(DDS, name=name, alpha=opt$pcutoff)
  #     res$mlog10p <- -log10(res$padj)
  #     write.table(data.frame(row_ID=rownames(res), res), file=file.path(opt$baseDir, opt$resultsDir, paste0(autoname, '_', name, '.lfc.tsv')),
  #                 sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  #   }
  #   
  #   # Shrunk LFC supposedly better for plotting and ranking
  #   for (name in coefficients[2:length(coefficients)]) {
  #     res <- lfcShrink(DDS, coef=name, type='apeglm')       # apeglm is the type recommended by DESeq2
  #     
  #     # Significance and counts are kept the same, but lfc is shrunk.
  #     res$mlog10p <- -log10(res$padj)
  #     write.table(data.frame(row_ID=rownames(res), res), file=file.path(opt$baseDir, opt$resultsDir, paste0(autoname, '_', name, '.shrunken.tsv')),
  #                 sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  #   }
  # }
  
}




