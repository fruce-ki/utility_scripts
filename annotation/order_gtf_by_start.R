library(data.table)

args <- commandArgs(trailingOnly=TRUE)
	
DT <- fread(args[1])

setorder(DT, V1, V4, V5)

fwrite(DT, file=args[1], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
