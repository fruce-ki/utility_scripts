#!/usr/bin/env Rscript

library(getopt)
library(data.table)

spec <- matrix(c(
  'controlGroups' , 'c', 1, "character", "Comma-separated list (no space) of group names for non-targeting sgRNAs.",
  'countsFile'    , 'f', 1, "character", "Tab-separated table of counts.",
  'groupCol'      , 'g', 2, "character", "Name of the column containing the group names ('group').",
  'help'          , 'h', 0, "logical"  , "Help.",
  'mincount'      , 'm', 2, "numeric"  , "Minimum count PER MILLION (will be scaled to libraries) to separate control guides into detected and non (1).",
  'guidesPerGene' , 'n', 2, "numeric"  , "Desired number of guides per control \"gene\" (6).",
  'outFile'       , 'o', 2, "character", "Output file (overwrite input file).",
  'nonrandom'     , 'r', 0, "logical"  , "Assign guides to groups in order of abundance instead of randomly (FALSE).",
  'seed'          , 's', 2, "numeric"  , "Seed for reproducible pseudo-randomisation.",
  'targetCol'     , 't', 2, "character", "Name of the column containing the guide names ('id').",
  'reference'     , 'z', 2, "character", "Comma separated column names across which to apply mincount for the controls. (if NULL, then all)"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)
# opt <- list(controlGroups='CTRLHS,CTRLMM', countsFile='/Volumes/groups/zuber/zubarchive/USERS/Kimon/markus/HLA1_staggered_v2/process/crispr-processed/counts/library/counts_mageck.txt', guidesPerGene=6, mincount=1, reference='NoIFNG_d0,NoIFNG_d0_2')

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

# Is/are there plasmid samples?
plasmids <- names(DT)[grepl('plasmid', names(DT))]

# For filtering, scale down the relevant samples to 1 million
DT2 <- as.data.table(lapply (DT[, union(plasmids, opt$reference), with=FALSE], function(x) { x / sum(x) * 1000000 }))

# Mark guides the are missing from plasmids or references
DT2[, plasmid_pass := rowSums(DT[, plasmids, with=FALSE] >= opt$mincount) >= 1] # if present in even just one plasmid
DT2[, ref_pass := rowSums(DT[, opt$reference, with=FALSE] >= opt$mincount) >= length(opt$reference)] # if present in all the references
message( paste(length(DT2$plasmid_pass) - sum(DT2$plasmid_pass), "guides removed for presence below", opt$mincount, "per million in", paste(plasmids, collapse=', '), ".") )

# Remove guides that are not present in the plasmids
DT <- DT[(DT2$plasmid_pass), ]
DT2 <- DT2[(plasmid_pass), ]

# Rename id and group columns to facilitate programming, try to avoid a name that might already exist.
names(DT)[which(names(DT)==opt$groupCol)] <- 'sexygroup666'
names(DT)[which(names(DT)==opt$targetCol)] <- 'sexyid666'

for (gr in groupNames){
  aux <- data.table(id = DT$sexyid666,
                    grp = DT$sexygroup666,
                    cnt = rowMeans(DT[, c(opt$reference), with=FALSE]),
                    nonzero = DT2$ref_pass,
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
