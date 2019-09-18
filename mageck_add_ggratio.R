#!/usr/bin/env Rscript

library(getopt)
library(data.table)
library(matrixStats)


spec = matrix(c(
  'infile'        , 'i', 1, "character", "Input gene results .tsv",
  'outfile'       , 'o', 1, "character", "Output .tsv",
  'gpat'          , 'p', 1, "character", "Pattern to recognise guides column names",
  'ggpat'         , 'q', 1, "character", "Pattern to recognise good_guides column names",
  'countsFile'    , 'f', 1, "character", "Tab-separated table of counts after nontargeting guides are grouped.",
  'mincount'      , 'm', 1, "numeric"  , "Minimum count to separate control guides into detected and non (1).",
  'reference'     , 'z', 1, "character", "Comma separated column names across which to apply mincount for the controls. (if NULL, then all)"
), byrow=TRUE, ncol=5)
opt = getopt(spec)
# opt <- list(infile='/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/HLA1_staggered_v2/process/mageck/genes_all_xref_dedup_reord_l10p_gg.tsv', countsFile='/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/HLA1_staggered_v2/process/crispr-processed/counts/library/counts_mageck_ctrls-grouped.txt', reference='NoIFNG_d0,NoIFNG_d0_2,plasmid', mincount=50)

if ( is.null(opt$gpat) )          { opt$gpat <- '\\.guides\\.' }
if ( is.null(opt$ggpat) )         { opt$ggpat <- '\\.guides_good\\.' }
if ( is.null(opt$targetCol) )     { opt$targetCol <- 'id' }
if ( is.null(opt$mincount) )      { opt$mincount <- 1 }
if (! is.null(opt$reference) )    { opt$reference <- unlist(strsplit(opt$reference, ',')) }


# Reference counts
DT <- fread(opt$countsFile, select=c("id", "group", opt$reference))

# Figure out the number of guides per group that pass the filtering threshold.
DT[ , pass := rowSums(DT[, -c(1,2)]) / length(opt$reference) >= opt$mincount ]
DT <- unique(DT[ , nguides := sum(pass), by=group] [, .(group, nguides)])
DT[nguides==0, nguides := NA_integer_]

# Input
DF <- fread(opt$infile, sep="\t", header=TRUE)

# Identify and isolate guides and good guides columns.
# selg <- grepl(opt$gpat, names(DF), perl=TRUE)
selgg <- grepl(opt$ggpat, names(DF), perl=TRUE)
selid <- grepl("group", names(DF), perl=TRUE)
# sel <- selgg | selg | selid
sel <- selgg | selid
subDF <- DF[, sel, with=FALSE]

# Calculate ratio
# # Rely on the order of columns out of the pipeline staying like this (likely).
# DF2 <- as.data.table(lapply(2:(length(subDF)/2 + 1), function(i) {      # going +1 past the last index to ensure the last pair is included, because the indices are backstepped from the loop value.
#   gg <- 2 * i -1
#   g <- gg -1
#   dif <- subDF[[gg]] / subDF[[g]]
#   return(dif)
# }))
subDF <- merge(subDF, DT, by="group", all.x=TRUE, all.y=FALSE)
DF2 <- subDF[, 2:(length(subDF)-1), with=FALSE] / subDF[, nguides]

DF2[is.na(DF2)] <- 0
# Rename
names(DF2) <- gsub(opt$ggpat, ".ggratio.", names(DF)[selgg], perl=TRUE)
# Add id for merging
DF2[, group := subDF$group]

fwrite(merge(merge(DF, DF2, by="group", all.x=TRUE), DT, by="group", all.x=TRUE),
       file=opt$outfile, sep="\t", col.names=TRUE, row.names=FALSE)
