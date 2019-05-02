#!/users/kimon.froussios/miniconda3/envs/mybasics/bin/Rscript

# Identify p-value columns as containing '.p.',
# compute their -log10(p)
# and output to a new file, indexed by 'id'.

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]  	# tsv table
outfile <- args[2] 	# output table for calc'ed values
ppat <- args[3]    	# pattern to recognise pvalue column names
id <- args[4]      	# name of index column

if (is.na(ppat))
	ppat <- '\\.p\\.'
if(is.na(id))
	id <- 'id'


library(data.table)

df <- fread(infile, sep="\t", header=TRUE)


# Identify and isolate id and pvalue-related columns.
sel <- (grepl(ppat, names(df), perl=TRUE) |
          grepl(paste0("^",id,"$"), names(df), perl=TRUE))
subdf <- df[, sel, with=FALSE]

# Transform
df2 = data.table(id = subdf[[id]])
for (n in names(subdf[, -c(id), with=FALSE])) {
  df2[[n]] = -log10(as.numeric(subdf[[n]]))
}
# Rename
names(df2) <- gsub(ppat, ".-Log10p.", names(df2), perl=TRUE)

fwrite(merge(df, df2, by=id, all.x=TRUE),
			file=outfile, sep="\t", col.names=TRUE, row.names=FALSE)
