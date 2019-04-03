#!/users/kimon.froussios/miniconda3/envs/mybasics/bin/Rscript

library(getopt)

spec = matrix(c(
  'controlGroups' , 'c', 1, "character", "Comma-separated list (no space) of group names for non-targeting sgRNAs.",
  'countsFile'    , 'f', 1, "character", "Tab-separated table of counts.",
  'groupCol'      , 'g', 2, "character", "Name of the column containing the group names ('group').",
  'help'          , 'h', 0, "logical"  , "Help.",
  'guidesPerGene' , 'n', 2, "numeric"  , "Desired number of guides per control \"gene\" (6).",
  'outFile'       , 'o', 2, "character", "Output file (overwrite input file).",
  'nonrandom'     , 'r', 0, "logical"  , "Assign guides to groups in order of appearance instead of randomly (FALSE).",
  'seed'          , 's', 2, "numeric"  , "Seed for reproducible pseudo-randomisation.",
  'targetCol'     , 't', 2, "character", "Name of the column containing the guide names ('id')."
), byrow=TRUE, ncol=5)
opt = getopt(spec)

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$guidesPerGene) )  { opt$guidesPerGene <- 6L }
if ( is.null(opt$outFile) )        { opt$outFile <- opt$countsFile }
if ( is.null(opt$targetCol) )      { opt$targetCol <- 'id' }
if ( is.null(opt$groupCol) )       { opt$groupCol <- 'group' }
if ( !is.null(opt$seed) )          { set.seed(opt$seed) }
if ( is.null(opt$nonrandom) )      { opt$nonrandom <- FALSE } else { opt$randomise <- TRUE }

# Group names
groupNames <- unlist(strsplit(opt$controlGroups, ',', fixed=TRUE))

library(data.table)

# Counts
dt <- fread(opt$countsFile)

# Rename id and group columns to facilitate programming, try to avoid a name that might already exist.
names(dt)[which(names(dt)==opt$groupCol)] <- 'sexygroup666'
names(dt)[which(names(dt)==opt$targetCol)] <- 'sexyid666'

for (gr in groupNames){
  # Group separately guides that are present andguides that are not present.
  nonzero <- rowSums(dt[, -c(1,2)]) > 0
  # Relevant rows
  sel <- dt$sexygroup666 == gr
  # Number of guides
  nx <- sum(sel & nonzero)
  nz <- sum(sel & !nonzero)
  # Number of subgroups to split into, decimals rounded up. The last group will have fewer guides if division was not whole.
  splitX <- floor(nx / opt$guidesPerGene) + 1
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
    newGroupsZ <- sample(newGroupsZ)
  }
  # And put the new values back into the table
  dt[sel & nonzero, sexygroup666 := newGroupsX]
  dt[sel & !nonzero, sexygroup666 := newGroupsZ]
}

# Restore original column names
names(dt)[which(names(dt)=='sexygroup666')] <- opt$groupCol
names(dt)[which(names(dt)=='sexyid666')] <- opt$targetCol

fwrite(dt, file=opt$outFile, sep = "\t", quote=FALSE)
