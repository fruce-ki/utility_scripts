#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec <- matrix(c(
  'controlGroups' , 'c', 1, "character", "Comma-separated list (no space) of group names for non-targeting sgRNAs.",
  'countsFile'    , 'f', 1, "character", "Tab-separated table of counts.",
  'groupCol'      , 'g', 1, "character", "Name of the column containing the group names ('group').",
  'help'          , 'h', 0, "logical"  , "Help.",
  'mincount'      , 'm', 1, "numeric"  , "Minimum count to separate control guides into detected and non (1).",
  'guidesPerGene' , 'n', 1, "numeric"  , "Desired number of guides per control \"gene\" (6).",
  'outFile'       , 'o', 1, "character", "Output file (overwrite input file).",
  'nonrandom'     , 'r', 0, "logical"  , "Assign guides to groups in order of abundance instead of randomly (FALSE).",
  'seed'          , 's', 1, "numeric"  , "Seed for reproducible pseudo-randomisation.",
  'targetCol'     , 't', 1, "character", "Name of the column containing the guide names ('id').",
  'reference'     , 'z', 1, "character", "Comma separated column names across which to apply mincount for the controls. (if NULL, then all)"#,
#  'ctrllist'      , 'l', 1, "character", "File where the list of control guide IDs will be written out for Mageck."
), byrow=TRUE, ncol=5)
opt <- getopt(spec)
# opt <- list(controlGroups='NONTARGETING', countsFile='/Volumes/groups/zuber/zubarchive/USERS/Kimon/sara_mfpl/TTP_gw_screen/process/preprocess/counts/library/counts_mageck.txt', guidesPerGene=6, mincount=50, reference='T0_Ref_1,T0_Ref_2')

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$guidesPerGene) ) { opt$guidesPerGene <- 6L }
if ( is.null(opt$outFile) )       { opt$outFile <- opt$countsFile }
if ( is.null(opt$targetCol) )     { opt$targetCol <- 'id' }
if ( is.null(opt$groupCol) )      { opt$groupCol <- 'group' }
if (! is.null(opt$seed) )         { set.seed(opt$seed) }
if ( is.null(opt$nonrandom) )     { opt$nonrandom <- FALSE } else { opt$randomise <- TRUE }
if ( is.null(opt$mincount) )      { opt$mincount <- 1 }
if (! is.null(opt$reference) )    { opt$reference <- unlist(strsplit(opt$reference, ',')) }

# Group 
groupNames <- unlist(strsplit(opt$controlGroups, ',', fixed=TRUE))

# Counts
DT <- fread(opt$countsFile)

# Rename id and group columns to facilitate programming, try to avoid a name that might already exist.
names(DT)[which(names(DT)==opt$groupCol)] <- 'sexygroup666'
names(DT)[which(names(DT)==opt$targetCol)] <- 'sexyid666'

for (gr in groupNames){
  # gr <- groupNames[1]
  
  # All guides above the threshold.
  meanCount <- NULL
  if (is.null(opt$reference)) {
    # All samples == columns excluding the guide and group ids.
    meanCount <- rowMeans(DT[, -c('sexygroup666', 'sexyid666'), with=FALSE])
  } else {
    # Specified samples only.
    meanCount <- rowMeans(DT[, c(opt$reference), with=FALSE])
  }
  
  aux <- data.table(id = DT$sexyid666,
                    grp = DT$sexygroup666,
                    cnt = meanCount,
                    nonzero = meanCount >= opt$mincount,
                    sel = DT$sexygroup666 == gr
          )[order(cnt)]
  aux <- aux[(sel), ]

  # Number of guides
  nx <- nrow(aux[(nonzero),])
  nz <- nrow(aux[(!nonzero),])
  # Shift a few guides around to ensure at least the
  # above-threshold group is divisible by the number of guides per gene.
  remainder <- nx %% opt$guidesPerGene
  if (remainder != 0){
    if ((nx - remainder > opt$guidesPerGene &&
         remainder <= opt$guidesPerGene / 2) ||
        nz < opt$guidesPerGene - remainder) {
      # If there are enough guides for at least one above-threshold group,
      # and the remainder is not big enough to justify rounding up,
      # or if there wouldn't be enough below-threshold guides to substract from,
      # reclassify the `remainder` bottom-count guides that are above threshold into below-threshold.
      aux[(nz + 1):(nz + remainder), nonzero := FALSE]
      nx <- nx - remainder
      nz <- nz + remainder

    } else {
      difference <- opt$guidesPerGene - remainder
      # Reclassify the `difference` top-count guides that are below threshold into above-threshold.
      aux[(nz - difference + 1):nz, nonzero := TRUE]
      # Round up.
      nx <- nx + difference
      nz <- nz - difference
    }
  }

  # Number of subgroups to split into, decimals rounded up. The last group will have fewer guides if division was not whole.
  splitX <- floor(nx / opt$guidesPerGene) + 1  # By now nx should be divisible, but allow freak case of nx < guidesPerGene.
  splitZ <- floor(nz / opt$guidesPerGene) + 1
  # Create new groupNames of suffixes
  newGroupsX <- paste0(gr, '_X_', rep(1:splitX, each = opt$guidesPerGene))
  newGroupsZ <- paste0(gr, '_Z_', rep(1:splitZ, each = opt$guidesPerGene))
  # Truncate extra spaces (due to having rounded up)
  newGroupsX <- newGroupsX[1:nx]
  newGroupsZ <- newGroupsZ[1:nz]
  # Shuffle for good measure
  if (!opt$nonrandom) {
    newGroupsX <- sample(newGroupsX)
    # newGroupsZ <- sample(newGroupsZ)
  }
  aux[(nonzero), grp := newGroupsX]
  aux[(!nonzero), grp := newGroupsZ]

  # And put the new values back into the table
  setkey(aux, id)
  setkey(DT, sexyid666)
  DT[aux$id, sexygroup666 := aux$grp]
}

# Restore original column names
names(DT)[which(names(DT)=='sexygroup666')] <- opt$groupCol
names(DT)[which(names(DT)=='sexyid666')] <- opt$targetCol

fwrite(DT, file=opt$outFile, sep = "\t", quote=FALSE)

# fileConn <- file('~/opt$outFile')
# writeLines(DT, fileConn)
# close(fileConn)
