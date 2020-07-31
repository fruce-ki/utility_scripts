#!/usr/bin/env Rscript

library(data.table)

# This script is not part of the pipeline. It is meant to be run stand-alone.
# The output is in .stats format that can be processed by mutpe

args <- commandArgs(trailingOnly = TRUE)
# args <- c('Ramos', '/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/merge_amplicons.txt')

vdj = args[1]
mergers <- fread(args[2], header=FALSE)

# Amplicons
if (vdj == 'B18') {
	amplicons <- data.frame(xleft = c(2022, 2195), xright=c(2209, 2381))
} else if (vdj == 'CH12') {
	amplicons <- data.frame(xleft = c(1694, 1874), xright=c(1874, 2077))
} else if (vdj == 'Ramos') {
	amplicons <- data.frame(xleft = c(3684, 3927, 4148), xright=c(3929, 4147, 4375))
} else {
	stop('Unknown view window')
}

# i <- 1
for (i in 1:nrow(mergers)) {
	files <- as.character(mergers[i,])
	outfile <- files[1]
	files <- files[2:length(files)]
	# files <- paste0('/Volumes', files)
	# Get rid of missing amplicon data
	files <- files[file.exists(files)]
	amplicons <- amplicons[file.exists(files), ]
	row.names(amplicons) <- files
	
	if (length(files) > 0){
		# Constrain each file to the amplicon coordinates (get rid of off-target alignments),
		# and then concatenate
		DT <- rbindlist(lapply(files, 
													 function(x) { 
													 		tmp <- fread(x)
													 		tmp[pos >= amplicons[x, "xleft"] & pos <= amplicons[x, "xright"]]
													 	}))
		setorder(DT, pos)
		# upgrade integer columns to float
		DT[, count := as.numeric(count)]
		DT[, depth := as.numeric(depth)]
		
		# Average overlapping positions (for the same nucleotide type)
		overlaps <- DT[duplicated(DT, by=c("pos", "type")), pos]
		DT[pos %in% overlaps,]
		DT[pos %in% overlaps, count := mean(count), by=c("pos", "type")]
		DT[pos %in% overlaps, depth := mean(depth), by=c("pos", "type")]
		# Lose the duplicate rows
		DT <- unique(DT)
	}
	
	fwrite(DT, file=outfile, quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}

