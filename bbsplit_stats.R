#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(getopt)

spec = matrix(c(
  'baseDir',      'b', 1, "character", "Full path to the project dir for nf-core/rnaseq.",
  'bbDir',        'd', 1, "character", "Relative sub-path to the bbsplit results directory."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

opt <- list(baseDir = '/SCRATCH/sBatch1/', bbDir = 'outputs/bbsplit/', sampleSheet = 'samplesheet.rnaseq.csv')

files <- dir(file.path(opt$baseDir, opt$bbDir))
files <- files[grepl('.stats.txt', files)]

DT <- Reduce(rbind,
             lapply(files, function(x) {
                            # x <- files[1]
                            tbl <- fread(file.path(opt$baseDir, opt$bbDir, x))
                            setnames(tbl, c("name", "unambiguousReadsPC", "unambiguousMB", "ambiguousReadsPC", "ambiguousMB",
                                           "unambiguousReads", "ambiguousReads", "assignedReads", "assignedBases"))

                            tbl[, sample := sub('.stats.txt', '', x)]
                            tbl
                          })
  )

DT[, row := 1:.N]

DT1 <- DT[, c("row", "sample", "name", "unambiguousReads", "ambiguousReads"), with=FALSE]
DT2 <- DT[, c("row", "sample", "name", "unambiguousReadsPC", "ambiguousReadsPC"), with=FALSE]

DT1 <- melt(DT1, id.vars = c("row", "sample", "name"), variable.name = "type", value.name = "PC", measure.vars = c("unambiguousReads", "ambiguousReads"))
DT2 <- melt(DT2, id.vars = c("row", "sample", "name"), variable.name = "type", value.name = "PC", measure.vars = c("unambiguousReadsPC", "ambiguousReadsPC"))

p1 <- ggplot(DT1, aes(x=sample, y=PC, fill=name)) +
  facet_wrap(~ type, nrow = 3, scales = "free_y") +
  geom_bar(stat="identity", position="stack") +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values=c("darkgreen", "gold3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

p2 <- ggplot(DT2, aes(x=sample, y=PC, fill=name)) +
  facet_wrap(~ type, nrow = 3, scales = "free_y") +
  geom_bar(stat="identity", position="stack") +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values=c("darkgreen", "gold3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


DT3 <- DT2[name == "primary" & type == "unambiguousReadsPC", ]
setkey(DT3, sample)

DT3[, bracket := 0]
DT3[PC >= 50, bracket := 50]
DT3[PC >= 30 & PC < 50, bracket := 30]
DT3[PC >= 20 & PC < 30, bracket := 20]
DT3[PC >= 15 & PC < 20, bracket := 15]
DT3[PC >= 5 & PC < 15, bracket := 5]
print('Brackets')
print(table(DT3$bracket))

DT3[, use := PC > 20]
print('Use')
print(table(DT3$use))

DT1 <- dcast(DT1[name == 'primary', ], sample ~ type, value.var = 'PC')
DT2 <- dcast(DT2[name == 'primary', ], sample ~ type, value.var = 'PC')

DT <- merge(DT1, DT2, by='sample')
DT <- merge(DT, DT3[, .(sample, bracket, use)], by = 'sample')

fwrite(DT, file = file.path(opt$baseDir, opt$bbDir, 'bbsplit.stats.csv'), sep=',', col.names = TRUE)

pdf(file = file.path(opt$baseDir, opt$bbDir, 'bbsplit.stats.pdf'), width = 10, height = 10)
print(p1)
print(p2)
dev.off()
