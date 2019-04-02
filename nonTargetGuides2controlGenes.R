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
  # Relevant rows
  sel <- dt$sexygroup666 == gr
  # Number of relevant guides
  n <- sum(sel)
  # Number of subgroups to split into, decimals rounded up. The last group will have fewer guides if division was not whole.
  splitN <- floor(n / opt$guidesPerGene) + 1
  # Create new groupNames of suffixes
  newGroups <- paste0(gr, '_', rep(1:splitN, each = opt$guidesPerGene))
  # Truncate extra spaces (due to having rounded up)
  newGroups <- newGroups[1:n]
  # Shuffle for good measure
  if (!opt$nonrandom)
    newGroups <- sample(newGroups)
  # And put the new values back into the table
  dt[sel, sexygroup666 := newGroups]
}

# Restore original column names
names(dt)[which(names(dt)=='sexygroup666')] <- opt$groupCol
names(dt)[which(names(dt)=='sexyid666')] <- opt$targetCol

fwrite(dt, file=opt$outFile, sep = "\t", quote=FALSE)
