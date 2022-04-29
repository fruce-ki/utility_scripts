library(data.table)

args <- commandArgs(trailingOnly=TRUE)
# args <- c("/Volumes/groups/busslinger/Kimon/tanja/R13065_RNAseq/process/duprates/duprates.txt")

DT <- fread(args[1])

DT1 <- melt(DT, id.vars="ID", measure.vars=names(DT)[grepl('dup[rR]ate', names(DT))], value.name="dupRate", variable.name="Sample")
DT1[, Sample := sub('.dup[rR]ate', '', Sample)]
setkey(DT1, ID, Sample)

DT2 <- melt(DT, id.vars="ID", measure.vars=names(DT)[grepl('RPKb', names(DT))], value.name="RPKb", variable.name="Sample")
DT2[, Sample := sub('.RPKb', '', Sample)]
setkey(DT2, ID, Sample)

DT <- merge(DT1, DT2)

fwrite(DT, file=sub('.txt$', '.melted.txt', args[1]), sep="\t", quote=FALSE)
