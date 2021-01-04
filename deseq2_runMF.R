#!/usr/bin/env Rscript

library(getopt)
library(DESeq2)


spec = matrix(c(
  'baseDir'      , 'b', 1, "character", "Base directory for everything (.)",
  'countsFile'   , 'f', 1, "character", "Tab-separated table of counts with `row_ID` and all the samples",
  'help'         , 'h', 0, "logical",   "Help",
  'resultsDir'   , 'o', 1, "character", "Directory in which to save the contrast results. If omitted, only the RDS will be output.",
  'RDSoutdir'    , 'r', 1, "character", "Directory in which to save the raw DESeq2 object (./process)",
  'samplesFile'  , 's', 1, "character", "Tab-separated table with `sample` column followed by the variable columns. Values must NOT contain punctuation other than '_' !!!",
  'control'      , 'x', 1, "character", "Comma separated value for each variable to use as reference for all comparisons, in the order variables are listed in samplesfile",
  'designFormula', 'd', 1, "character", "Design formula",
  'reducedFormula','R', 1, "character", "Reduced formula for LR test",
  'forvar',        'F', 1, "character", "Looping variable (ie. carry out the design formula for each level of this var, separately).",
  'selvar',        'S', 1, "character", "Selection variable (always together with sellev). It is applied BEFORE forvar.",
  'sellev',        'L', 1, "character", "Selection variable level. Together with selvar, it is used to reduce the samplesFile to just the rows having that value in that variable.",
  'minCount'     , 'M', 1, "numeric",   "Minimum number of reads across all samples combined for a gene to be considered (10).",
  'reportTemplate','T', 1, "character", "Template Rmd file for DE report.",
  'pcutoff',       'p', 1, "numeric",   "P-value cutoff (0.05)",
  'ntop',          'n', 1, "numeric",   "Number of hits to highlight (50)",
  'lfcthreshold',  'c', 1, "numeric",   "Log2 fold-change threshold (1)",
  'nidcols',       'I', 1, "numeric",   "Number of ID columns at the start of the table (1). see also -i",
  'idcol',         'i', 1, "numeric",   "ID column to use (1). The others will be removed. See also -I.",
  'bmF',           'B', 0, "logical",   "Disable baseMean filtering. (False)",
  'all',           'A', 0, "logical",   "Enable round-robin-style pairwise comparisons. Otherwise only comparisons only against reference. (False)" 
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(baseDir='/Volumes/groups/zuber/zubarchive/USERS/Kimon/jakub/M10716_slamseq', countsFile='process_nonmerged/slam/all_collapsed_readCount-only_xref.txt', resultsDir='results_nonmerged/DE', RDSoutdir='process_nonmerged/DE', samplesFile='description/covars_nonmerged.txt', control='THP1,WT,DMSO,HNWTTDRXX,HNWTTDRXX_2', designFormula='~ treatment', forvar='host', minCount=10, lfcthreshold=1, nidcols=3, idcol=2, ntop=50, bmF=FALSE, pcutoff=0.05, bmF=FALSE, all=FALSE)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$reportFile) ) {
  opt$reportFile <- '~/utility_scripts/deseq2_report_template.Rmd'
}

if ( is.null(opt$baseDir    ) ) { 
  opt$baseDir <- '.'         
}
if ( is.null(opt$RDSoutdir ) ) { 
  opt$RDSoutdir <- './process' 
}

if ( is.null(opt$minCount ) ) { 
  opt$minCount <- 10
}
if ( is.null(opt$ntop ) ) { 
  opt$ntop <- 50
}

if ( is.null(opt$pcutoff ) ) { 
  opt$pcutoff <- 0.05
}
if ( is.null(opt$lfcthreshold ) ) { 
  opt$lfcthreshold <- 1
}

if ( is.null(opt$nidcols ) ) { 
  opt$nidcols <- 1
}
if ( is.null(opt$idcol ) ) { 
  opt$idcol <- 1
}

if ( is.null(opt$bmF ) ) { 
  opt$bmF <- FALSE
}
if ( is.null(opt$all ) ) { 
  opt$all <- FALSE
}
opt$bmF <- !opt$bmF    # Convert -B flag functionality from switch-on to switch-off. 



# Output destinations
dir.create(file.path(opt$baseDir, opt$RDSoutdir), recursive=TRUE)

if ( !is.null(opt$resultsDir ) ) { 
  dir.create(file.path(opt$baseDir, opt$resultsDir), recursive=TRUE)
}


# Input
covars <- as.data.frame(read.csv(file.path(opt$baseDir, opt$samplesFile), sep="\t", header=TRUE, row.names='sample', check.names = FALSE))
cts <- read.csv(file.path(opt$baseDir, opt$countsFile), sep="\t", header=TRUE, check.names = FALSE) 
# cts <- cts[, c(names(counts)[1:opt$nidcols], covars$sample)]    # Only samples specified in covars. ## It is handled somewhere in the RMD. If done here it causes a crash!

# Cut out extra annotation/id columns
xref <- cts[, 1:opt$nidcols]
cts <- round(as.matrix(cts[, (opt$nidcols+1):ncol(cts)]), digits=0)
rownames(cts) <- xref[, opt$idcol]

# R adds 'X' to the begininng fo matric colnames that start with a digit, regardless of whether they can be misinterpreted as a number or not.
# That breaks the association with the names in covars and DESeq2 crashes. So the X needs to be added in covars too.

# Calculate RPMs. 
# (for 3'seq like quantseq and slamseq, length normalisation is not approriate)
colsums <- colSums(cts, na.rm=TRUE)
RPMs <- as.matrix(as.data.frame(lapply(colnames(cts), function(n) { cts[,n] / colsums[n] * 1e6 }) ))
colnames(RPMs) <- colnames(cts)
rownames(RPMs) <- rownames(cts)

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
  # in case not all combinatinons of variables exist
  covars[[opt$forvar]] <- droplevels(covars[[opt$forvar]])
  
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
  subrpm <- RPMs[, match(cdn, colnames(RPMs))]
  subrpm[is.na(subrpm)] <- 0

  # DESeq2 begins here
  DDS <- DESeqDataSetFromMatrix(countData = subcts,
                                colData = coldata,
                                design = as.formula(opt$designFormula))
  # Set reference condition(s)
  for (n in vars){
    DDS[[n]] <- relevel(DDS[[n]], ref=refs[n])
  }
  
  # Prefilter out genes with very poor coverage.
  # This is a speed and memory improvement. DESeq2 later internally applies additional thresholds.
  keep <- rowSums(counts(DDS)) >= opt$minCount
  DDS <- DDS[keep, ]
  subrpm <- subrpm[keep, ]
  # subxref <- xref[keep, ]
  
  # Prepare file name.
  autoname <- NULL
  if (!is.null(opt$selvar)) {
    autoname <- paste(opt$selvar, opt$sellev, sep="-")
  }
  if (!is.null(opt$forvar)) {
    autoname <- paste(autoname, paste(opt$forvar, V, sep="-"), sep="_")
  }
  
  ddsrds <- file.path(opt$baseDir, opt$RDSoutdir, paste0(paste(autoname, gsub(' |~', '', opt$designFormula), sep='_'), '_deseq2data.RDS'))
  saveRDS(DDS, file=ddsrds)
  rpmrds <- file.path(opt$baseDir, opt$RDSoutdir, paste0(paste(autoname, gsub(' |~', '', opt$designFormula), sep='_'), '_rpm.RDS'))
  saveRDS(subrpm, file=rpmrds)
  
  # Fire up an Rmd report
  rmarkdown::render(opt$reportFile, 
                    # output_file = paste0(paste(autoname, name, sep='_'), '_deseq2.html'), 
                    output_file = paste0(paste(autoname, gsub(' |~', '', opt$designFormula), sep='_'), '_deseq2.html'),
                    output_dir = file.path(opt$baseDir, opt$resultsDir),
                    params=list(ntop = opt$ntop,
                                pcutoff = opt$pcutoff,
                                lfcthreshold = opt$lfcthreshold,
                                baseDir = opt$baseDir,
                                resultsDir = opt$resultsDir,
                                prefix = autoname,
                                derds = ddsrds,
                                rmrds = rpmrds,
                                minCount = opt$minCount,
                                filterBaseMean = opt$bmF,
                                reducedFormula = opt$reducedFormula,
                                roundrobin=opt$all)
                    )
  unlink(rpmrds)
  
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




