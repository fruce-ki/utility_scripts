#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/zuber/zubarchive/USERS/Kimon/sara_mfpl/TTP_gw_screen/process/mageck/guides_all_xref.tsv', 'group')

file <- args[1]
colname <- args[2]

DT <- fread(file)
dups <- which(grepl(colname, names(DT), fixed=TRUE))

merged <- apply(DT[, grepl(colname, names(DT), fixed=TRUE), with=FALSE], 1, function(x) {
	y=unique(x)
	if (all(is.na(y))) {
		NA_character_
	} else {
		y[!is.na(y)][1]
	}
})

DT <- DT[, ! 1:length(DT) %in% dups[2:length(dups)], with=FALSE]
set(DT, i=NULL, j=colname, merged)

fwrite(DT, file=sub('\\.tsv$|\\.txt$', '_dedup.tsv', file, perl=TRUE), sep="\t", quote=FALSE)
