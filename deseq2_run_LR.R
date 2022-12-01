#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'bmF',           'B', 0, "logical",   "Disable baseMean filtering. (False)",
  'baseDir',       'b', 1, "character", "Base directory for everything (.)",
  'lfcthreshold',  'c', 1, "numeric",   "Log2 fold-change threshold (1)",
  'fc4p',          'C', 0, "logical",   "Incorporate the LFC threshold in the Null Hypothesis.",
  'createID',      'D', 0, "logical",   "Non-summarized intron or exon counts by featureCounts. The first column is not unique IDs because genes are repeated.",
  'countsFile',    'f', 1, "character", "Tab-separated table of raw counts with `row_ID` and all the samples.",
  'fullFormula',   'F', 1, "character", "Full model formula for LR or Wald test",
  'specnorm',      'G', 1, "character", "Special normalisation. Gene ids that match this pattern will not be included for calculation of library sizes for TPM/RPM normalisation.",
  'help',          'h', 0, "logical",   "This help.",
  'altHypo',       'H', 1, "character", "Alternative hypothesis type (greaterAbs).",
  'idcol',         'i', 1, "numeric",   "ID column to use (1). The others will be removed. See also -I.",
  'nidcols',       'I', 1, "numeric",   "Number of ID columns at the start of the table (1). See also -i.",
  'prescaled',     'k', 0, "logical",   "Don't let Deseq2 rescale the libraries. (False)",
  'label',         'l', 0, "logical",   "Add comparison details to the standard column names. Useful for merging multiple outputs.",
  'minCount',      'm', 1, "numeric",   "Minimum number of reads in either condition in a contrast. (100)",
  'minTPM',        'M', 1, "numeric",   "Minimum TPM in either condition in a contrast. (5)",
  'ntop',          'n', 1, "numeric",   "Number of hits to highlight (50)",
  'resultsDir',    'o', 1, "character", "Directory in which to save the contrast results. If omitted, only the RDS will be output.",
  'pcutoff',       'p', 1, "numeric",   "P-value cutoff (0.05)",
  'prefix',        'P', 1, "character", "Prefix to add to output file name and column names (''). ",
  'RDSoutdir',     'r', 1, "character", "Directory in which to save the raw DESeq2 object (./process)",
  'reducedFormula','R', 1, "character", "Reduced formula for LR test.",
  'samplesFile',   's', 1, "character", "Tab-separated table with `Sample` column followed by the variable columns. Values must NOT contain punctuation other than '_' !!!",
  'reportTemplate','T', 1, "character", "Template Rmd file for DE report.",
  'vst',           'v', 0, "logical",   "Plot VST counts instead of rlog2",
  'widthsFile',    'w', 1, "character", "Full path to table of feature lengths",
  'widthsCol',     'W', 1, "integer",   "Column in the counts file that contains feature lengths. This should be among the non-count columns, see also -I. It is an alternative to providing a separate lengths file.",
  'comparisons',   'x', 1, "character", "Comma separated list indicating which levels of the formula factor to compare, in order of appearance. ie 'Treatment-2v1,Treatment-3v1,Treatment-3v2,Batch-4v5'.",
  'contexts',      'X', 1, "character", "Comma separated list of covariate-level pairs indicating the context/subset to use for each comparison, ie 'celltype-A,celltype-A,celltype-B,culture-Z"
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(createID=FALSE, baseDir='/scratch-cbe/users/kimon.froussios/tanja/R14425_RNAseq', countsFile='process/featureCounts_fixed/exon_spliced_genecounts.txt', resultsDir='results/DE', samplesFile='description/covars.txt', minCount=0, minTPM=0, lfcthreshold=1, nidcols=6, idcol=1, ntop=50, bmF=FALSE, pcutoff=0.05, prescaled=FALSE, prefix='exon_spliced_genecounts', label=FALSE, widthsCol=6, reportTemplate="/groups/busslinger/Kimon/tanja/R14425_RNAseq/code/deseq2_report_template_LR.Rmd", comparisons="TF_type-1v2,TF_type-3v4,TF_type-3v4,TF_type-3v4,TF_type-3v4,TF_type-3v4,TF_type-5v6", contexts="Cell-1,Cell-1,Cell-2,Cell-3,Cell-4,Cell-5,Cell-5", fullFormula="~TF_type", reducedFormula="NULL", specnorm='^ENSMUSG|^TCR|^IGHG|^IGHM|^IGH[^GM]')


if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}



library(DESeq2)
library(data.table)



stopifnot(!is.null(opt$countsFile))
stopifnot(!is.null(opt$samplesFile))

if (is.null(opt$reportTemplate))  opt$reportTemplate <- '~/utility_scripts/deseq2_report_template_SF.Rmd'
if (is.null(opt$createID))  opt$createID <- FALSE 
if (!is.null(opt$specnorm) && opt$specnorm == "NULL") opt$specnorm <- NULL
if (is.null(opt$prescaled)) opt$prescaled <- FALSE
if (is.null(opt$label)) opt$label <- FALSE
if (is.null(opt$baseDir)) opt$baseDir <- '.'
if (is.null(opt$resultsDir))  opt$resultsDir <- './results' 
if (is.null(opt$minCount)) opt$minCount <- 100
if (is.null(opt$minTPM)) opt$minTPM <- 5
if (is.null(opt$ntop))  opt$ntop <- 50
if (is.null(opt$pcutoff)) opt$pcutoff <- 0.05
if (is.null(opt$lfcthreshold))  opt$lfcthreshold <- 1
if (is.null(opt$altHypo)) opt$altHypo <- 'greaterAbs'
if (is.null(opt$nidcols)) opt$nidcols <- 1
if (is.null(opt$idcol)) opt$idcol <- 1
if (is.null(opt$bmF)) opt$bmF <- FALSE
opt$bmF <- !opt$bmF    # Convert -B flag functionality from switch-on to switch-off. 
if (is.null(opt$vst)) opt$vst <- FALSE
if (!is.null(opt$fullFormula) && opt$fullFormula == "NULL") opt$fullFormula <- NULL
if (!is.null(opt$reducedFormula) && opt$reducedFormula == "NULL") opt$reducedFormula <- NULL
if (!is.null(opt$comparisons) && opt$comparisons == "NULL") opt$comparisons <- NULL
if (is.null(opt$contexts) | opt$contexts == "NULL") stop("Need at least one context!")

# Output destinations
dir.create(file.path(opt$baseDir, opt$resultsDir, opt$prefix), recursive=TRUE)

# Covariates
covars <- read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", row.names="Sample", header=TRUE, check.names = FALSE, colClasses="character")

# Contrasts
contexts <- strsplit(gsub(' ', '', opt$contexts), ',')[[1]]
contexts <- strsplit(contexts , '-')
if (!is.null(opt$comparisons)) {
  comparisons <- strsplit(gsub(' ', '', opt$comparisons), ',')[[1]]
  
  stopifnot(length(contexts) == length(comparisons))
  
  comparisons <- strsplit(strsplit(opt$comparisons, ',')[[1]], '-')
  comparisons <- lapply(comparisons, function(x){
    # x <- comparisons[[1]]
    list( x[1], as.integer(strsplit(x[2], 'v')[[1]]) )
  })
}

# Counts
cts <- fread(file.path(opt$baseDir, opt$countsFile), sep="\t", header=TRUE, check.names = FALSE, colClasses="character")
    # Setting all columns to character prevents numeric IDs at the top of the ID columns from auto-defining these columns as numeric when there could be character IDs further down that would otherwise become NaN.
cts <- unique(cts)         # non-summarized featureCounts have repeated exons from sharing among transcripts.
# Remove disruptive features (such as immunoglobins in plasmablasts).
ogn <- nrow(cts)
if (!is.null(opt$specnorm))
  cts <- cts[!grepl(opt$specnorm, Geneid), ]
deltan <- ogn - nrow(cts)
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


#######################
#######################
###  Split up data as instructed
for(y in unique(contexts)) {
  # y <- unique(contexts)[[3]]
  
  contvar <- y[1]
  contlev <- unique(covars[[contvar]])[as.integer(y[2])]
  
  sel <- which(vapply(contexts, function(a, b=y) { a[1] == b[1] && a[2] == b[2] },  logical(1)))
  
  condvars <- NULL
  treatlev <- NULL
  reflevs <- NULL
  if (!is.null(opt$comparisons)) {
    selcomps <- comparisons[sel]
    condvars <- vapply(selcomps, function(x){ x[[1]] }, character(1)) 
    
    stopifnot(length(unique(condvars)) == 1)
    
    treatlevs <- vapply(selcomps, function(x){ as.character(unique(covars[[x[[1]]]])[x[[2]][1]]) }, character(1))
    reflevs <- vapply(selcomps, function(x){ as.character(unique(covars[[x[[1]]]])[x[[2]][2]]) }, character(1))
    
    sanity <- vapply(condvars, function(x) { grepl(x, opt$fullFormula, fixed=TRUE) }, logical(1))
    if (!is.null(opt$fullFormula) && any(!sanity)) {
      warning(paste(condvars[which(!sanity)], "not in formula", opt$fullFormula))
    }
  }
  
  # Remove unused samples from the counts table,
  # and ensure samples are in the same order as the covariates table.
  subcovars <- covars[ covars[[contvar]] == contlev, ]
  
  if (!is.null(opt$comparisons)) {
    stopifnot(all( vapply(1:length(condvars), function(i) { reflevs[i] %in% subcovars[[condvars[i]]] }, logical(1)) ))
    stopifnot(all( vapply(1:length(condvars), function(i) { treatlevs[i] %in% subcovars[[condvars[i]]] }, logical(1)) ))
  }
  
  cdn <- rownames(subcovars)
  ctsn <- colnames(cts)
  subcts <- cts[, ctsn[match(cdn, ctsn)]]
  # Fill in non-detections
  subcts[is.na(subcts)] <- 0

  # DESeq2 begins here
  designFormula <- ifelse(!is.null(opt$fullFormula),
                                   opt$fullFormula,
                                   paste0('~', paste(unique(condvars),collapse='+')) )
  DDS <- DESeqDataSetFromMatrix(countData = round(subcts, digits=0),
                                colData = subcovars,
                                design = as.formula(designFormula) )
                                
  # To scale or not to scale
  if (opt$prescaled) {
    factors <- rep(1, times=length(cdn))
    names(factors) <- cdn
    sizeFactors(DDS) <- factors
  }

  
  # Prepare file name.
  if (opt$prefix == '') {
    # prefix <- paste(paste(contvar, contlev, sep='_'), designFormula, sep='.')
    prefix <- paste(contvar, contlev, sep='_')
  } else {
    # prefix <- paste(opt$prefix, paste(contvar, contlev, sep='_'), designFormula, sep='.')
    prefix <- paste(opt$prefix, paste(contvar, contlev, sep='_'), sep='.')
  }
  name <- paste(prefix, 'deseq2.html', sep='.')
  
  #######################
  #######################
  
  subcomps <- NULL
  if (!is.null(opt$comparisons)) {
    subcomps <- lapply(1:length(condvars), function(i){ c(condvars[i], treatlevs[i], reflevs[i]) })
  }
  
  # Fire up an Rmd report
  rmarkdown::render(opt$reportTemplate,
                    output_file = name,
                    output_dir = file.path(opt$baseDir, opt$resultsDir, opt$prefix),
                    intermediates_dir=file.path(opt$baseDir, opt$resultsDir, sub('deseq2.html', 'tmp', name)),
                    params=list(selfname = name,
                                ntop = opt$ntop,
                                pcutoff = opt$pcutoff,
                                lfcthreshold = opt$lfcthreshold,
                                baseDir = opt$baseDir,
                                resultsDir = file.path(opt$resultsDir, opt$prefix),
                                prefix = prefix,
                                dds = DDS,
                                minCount = opt$minCount,
                                minTPM = opt$minTPM,
                                filterBaseMean = opt$bmF,
                                reducedFormula = opt$reducedFormula,
                                longlabel=opt$label,
                                covars=covars,
                                widths = featLens,
                                specnorm=opt$specnorm,
                                deltaN=deltan,
                                rlog=!opt$vst,
                                comparisons=subcomps,
                                althyp=opt$altHypo,
                                fc4p=opt$fc4p)
  )
}



