#!/usr/bin/env Rscript

library(getopt)
library(dupRadar)
library(data.table)

spec = matrix(c(
	'bam',            'b', 1, "character", "Full path to a BAM file. Must have duplicate reads already marked by picard or bamutils.",
	'gtf',            'g', 1, "character", "Full path to a GTF file.",
	'output',         'o', 1, "character", "Full path to a destination file.",
	'paired',         'p', 0, "logical",   "Paired-end sequencing (FALSE).",
	'stranded',       's', 1, "integer",   "Stranded sequencing. 0 unstranded, 1 stranded, 2 reversely stranded (0).",
	'threads',        't', 1, "integer",   "Threads (1)."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(bam="/scratch-cbe/users/kimon.froussios/robyn/R12658_RNAseq/data_chrfixed/Tcf3_noTir1_Auxin_3_5_m__180007.aligned.bam", gtf="/groups/busslinger/Kimon/genome_mm10/refGene.mm10.2018_1206.withIgtcr.final.annotation.IMP_IGTCR_only.with_introns.gtf", paired=TRUE, stranded=1, threads=1, minCount=20, minRPKM=5, ntop=25, output='/scratch-cbe/users/kimon.froussios/robyn/R12658_RNAseq/duprate/Tcf3_noTir1_Auxin_3_5_m__180007.duprate.txt')

if (is.null(opt$bam)) stop('No BAM input')
if (is.null(opt$gtf)) stop('No GTF input')
if (is.null(opt$output)) stop('No output')
if (is.null(opt$stranded)) opt$stranded <- 0
if (is.null(opt$paired)) opt$paired <- FALSE
if (is.null(opt$threads)) opt$threads <- 1


wd <- setwd(dirname(opt$output))


dm <- as.data.table(analyzeDuprates(opt$bam, opt$gtf, opt$stranded, opt$paired, opt$threads))
setorder(dm, -dupRate, na.last = TRUE)
fwrite(dm, file=opt$output, quote=FALSE, sep="\t")


subdm <- dm[filteredCounts >= opt$minCount & RPKM >= opt$minRPKM, ]
# subdm <- dm[RPKM >= opt$minRPKM, ]
setorder(subdm, -dupRate, na.last = TRUE)


pdf(sub('txt$|tsv$', 'pdf', opt$output))

plot.new()
text(x=.1, y=.1, basename(opt$bam))

print( duprateExpDensPlot(DupMat=dm) )
print( duprateExpBoxplot(DupMat=dm) )
print( readcountExpBoxplot(DupMat=dm) )

dev.off()


setwd(wd)
