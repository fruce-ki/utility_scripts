#!/usr/bin/env Rscript

# Identify p-value columns as containing '.p.',
# compute their -log10(p)
# and output to a new file, indexed by 'id'.

library(getopt)

spec = matrix(c(
  'infile'        , 'i', 1, "character", "Input .tsv",
  'outfile'       , 'o', 1, "character", "Output .tsv",
  'ppat'          , 'p', 1, "character", "Pattern to recognise pvalue column names",
  'id'            , 'r', 1, "character", "Column name for row IDs"
), byrow=TRUE, ncol=5)
opt = getopt(spec)

if (is.null(opt$ppat))
  opt$ppat <- '\\.p\\.'
if(is.null(opt$id))
  opt$id <- 'id'


library(data.table)

df <- fread(opt$infile, sep="\t", header=TRUE)


# Identify and isolate id and pvalue-related columns.
sel <- (grepl(opt$ppat, names(df), perl=TRUE) |
          grepl(paste0("^",opt$id,"$"), names(df), perl=TRUE))
subdf <- df[, sel, with=FALSE]

# Transform
df2 <- data.table(id = subdf[[opt$id]])
for (n in names(subdf[, -c(opt$id), with=FALSE])) {
  df2[[n]] = -log10(as.numeric(subdf[[n]]))
}
# Rename
names(df2) <- gsub(opt$ppat, ".-Log10p.", names(df2), perl=TRUE)

fwrite(merge(df, df2, by.x=opt$id, by.y="id", all.x=TRUE),
			file=opt$outfile, sep="\t", col.names=TRUE, row.names=FALSE)
