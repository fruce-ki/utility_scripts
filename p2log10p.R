#!/users/kimon.froussios/miniconda3/envs/mybasics/bin/Rscript

# Identify p-value columns as containing '.p.', 
# compute their -log10(p)
# and output to a new file, indexed by 'id'.

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
ppat <- args[3]
id <- args[4]

if (is.na(ppat))
	ppat <- '\.p\.'
if(is.na(id))
	id <- 'id'

library(data.table)

df <- fread(infile, sep="\t", header=TRUE)

# Identify and isolate id and pvalue-related columns.
sel <- (grepl(ppat, names(df), fixed=FALSE) | 
          grepl(paste0("^",id,"$"), names(df), fixed=FALSE)) 
df <- df[, sel, with=FALSE]

# Transform
df2 = data.table(id = df[[id]])
for (n in names(df[, -c(id)])) {
  df2[[n]] = -log10(as.numeric(df[[n]]))
}
# Rename
names(df2) <- gsub(ppat, ".-Log10p.", names(df2), fixed=FALSE)

fwrite(df2, file=outfile, sep="\t", col.names=TRUE, row.names=FALSE)

