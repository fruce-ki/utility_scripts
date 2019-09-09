#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'infile'        , 'i', 1, "character", "Input .tsv",
  'outfile'       , 'o', 1, "character", "Output .tsv",
  'gpat'          , 'p', 1, "character", "Pattern to recognise guides column names",
  'ggpat'         , 'q', 1, "character", "Pattern to recognise good_guides column names",
  'id'            , 'r', 1, "character", "Column name for row IDs"
), byrow=TRUE, ncol=5)
opt = getopt(spec)

if (is.null(opt$gpat))
  opt$gpat <- '\\.guides\\.'
if (is.null(opt$ggpat))
  opt$ggpat <- '\\.guides_good\\.'

if(is.null(opt$id))
  opt$id <- 'id'


library(data.table)

df <- fread(opt$infile, sep="\t", header=TRUE)

# Identify and isolate guides and good guides columns.
selg <- grepl(opt$gpat, names(df), perl=TRUE)
selgg <- grepl(opt$ggpat, names(df), perl=TRUE)
selid <- grepl(paste0("^",opt$id,"$"), names(df), perl=TRUE)
sel <- selgg | selg | selid
subdf <- df[, sel, with=FALSE]

# Calculate ratio
# Rely on the order of columns out of the pipeline staying like this (likely).
df2 <- as.data.table(lapply(2:(length(subdf)/2 + 1), function(i) {      # going +1 past the last index to ensure the last pair is included, because the indices are backstepped from the loop value.
  gg <- 2 * i -1
  g <- gg -1
  dif <- subdf[[gg]] / subdf[[g]]
  return(dif)
}))
df2[is.na(df2)] <- 0
# Rename
names(df2) <- gsub(opt$ggpat, ".ggratio.", names(df)[selgg], perl=TRUE)
# Add id for merging
df2[, id := subdf[[opt$id]]]

fwrite(merge(df, df2, by.x=opt$id, by.y="id", all.x=TRUE),
       file=opt$outfile, sep="\t", col.names=TRUE, row.names=FALSE)
