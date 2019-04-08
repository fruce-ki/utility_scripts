#!~/miniconda3/envs/mybasics/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
cols <- unlist(strplit(args[3], ','))   # comma separated list of column names (or numeric indices if header-less)
newname <- args[4]                      # if headerless input, this colname won't show in the output
replace <- as.logical(args[5])
header <- as.logical(args[6])

library(data.table)

df <- fread(infile, header=header)
if (! header)
  cols <- as.nu
df[, newname, with=FALSE] <- rowSums(df[, cols, with=FALSE])
if (replace) {
  df[, cols, with=FALSE] <- NULL
}
fwrite(df, file=outfile, header=header)
