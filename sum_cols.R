#!~/miniconda3/envs/mybasics/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
cols <- unlist(strsplit(args[3], ','))   # comma separated list of column names (or numeric indices if header-less)
newname <- args[4]                      # if headerless input, this colname won't show in the output
replace <- as.logical(args[5])
header <- as.logical(args[6])

library(data.table)

df <- fread(infile, header=header)

# In a headerless file, cols must be positional indices.
# In a headered file, cols must be given by name.
if (! header)
  cols <- as.numeric(cols)

# Sum. Give a temporary new name to the new column, to prevent clash with existing columns in case of replace==TRUE.
v <- rowSums(df[, cols, with=FALSE])
df[, newsexycolumnname666 := v]

# Drop the old columns.
if (replace) {
  df[, cols] <- NULL
}
# Now name the new column properly.
names(df)[length(names(df))] <- newname

fwrite(df, file=outfile, col.names=header, row.names=FALSE)
