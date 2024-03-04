#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(patchwork)
library(getopt)

spec = matrix(c(
  'baseDir',      'b', 1, "character", "Full path to the project dir for nf-core/rnaseq.",
  'bbDir',        'd', 1, "character", "Relative sub-path to the bbsplit results directory.",
  'pc_thresh',    'p', 1, "numeric", "Minimum acceptable percentage of human content (20).",
  'cnt_thresh',   'c', 1, "numeric", "Minimum acceptable read count for human content (5000000)."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

# opt <- list(baseDir = '/SCRATCH/sBatch1/', bbDir = 'outputs/bbsplit/', sampleSheet = 'samplesheet.rnaseq.csv')

files <- dir(file.path(opt$baseDir, opt$bbDir))
files <- files[grepl('.stats.txt', files)]

if (is.null(opt$pc_thresh)) opt$pc_thresh <- 20
if (is.null(opt$cnt_thresh)) opt$cnt_thresh <- 5e6

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

DT1 <- melt(DT1, id.vars = c("row", "sample", "name"), variable.name = "type", value.name = "Count", measure.vars = c("unambiguousReads", "ambiguousReads"))
DT2 <- melt(DT2, id.vars = c("row", "sample", "name"), variable.name = "type", value.name = "PC", measure.vars = c("unambiguousReadsPC", "ambiguousReadsPC"))

# Ambiguous are counted doubly
DT1 <- DT1[!(name != 'primary' & type == 'ambiguousReads'), ]
DT2 <- DT2[!(name != 'primary' & type == 'ambiguousReadsPC'), ]
DT1[type == 'ambiguousReads', name := 'ambiguous']
DT2[type == 'ambiguousReadsPC', name := 'ambiguous']

DT1[, name := factor(name, ordered = TRUE, levels = rev(c('primary', sort(unique(DT1[name != 'primary', name])))))]
DT2[, name := factor(name, ordered = TRUE, levels = rev(c('primary', sort(unique(DT2[name != 'primary', name])))))]


fills <- c('orange3', 'darkgreen', rep('yellow3', times = length(levels(DT1$name)) - 2))
names(fills) <- c('ambiguous', 'primary', levels(DT1$name)[1:(length(levels(DT1$name)) - 2)])


p1 <- ggplot(DT1, aes(x=sample, y=Count, fill=name)) +
  geom_bar(stat="identity", position="stack") +
  geom_hline(, yintercept = opt$cnt_thresh, linetype = 'dashed', colour = 'blue') +
  scale_x_discrete(position = "bottom") +
  scale_fill_manual(values=fills) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  #theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

p2 <- ggplot(DT2, aes(x=sample, y=PC, fill=name)) +
  geom_bar(stat="identity", position="stack") +
  geom_hline(yintercept = opt$pc_thresh, linetype = 'dashed', colour = 'blue') +
  scale_x_discrete(position = "top") +
  #scale_y_reverse() +
  scale_fill_manual(values=fills) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, vjust=0.5), axis.title.x = element_blank())




DT3 <- merge(DT1[name == "primary" & type == "unambiguousReads", ],
             DT2[name == "primary" & type == "unambiguousReadsPC", ],
             by = c('sample', 'name', 'row'))
setnames(DT3, c('sample', 'drop0', 'drop1', 'drop2', 'unambiguousReads', 'drop3', 'unambiguousReadsPC'))
DT3 <- DT3[, !grepl('drop', names(DT3)), with = FALSE]
setkey(DT3, sample)



brackets <- sort(unique(c(opt$cnt_thresh, 1e7, 50e6, 20e6, 10e6, 5e6, 2e6, 0)), decreasing = FALSE)  # include large enough top limit and 0 for the conditional to work correctly
for (i in 1:(length(brackets) - 1) ) {
  DT3[unambiguousReads >= brackets[i] & unambiguousReads < brackets[i+1], bracket := paste0('<', brackets[i+1] / 1e6, 'M')]
}


bracketsPC <- sort(unique(c(opt$pc_thresh, 100, 75, 50, 35, 20, 10, 0)), decreasing = FALSE)  # include 100 and 0 for the conditional to work correctly
for (i in 1:(length(bracketsPC) - 1) ) {
  DT3[unambiguousReadsPC >= bracketsPC[i] & unambiguousReadsPC < bracketsPC[i+1], bracketPC := paste0('<', bracketsPC[i+1], '%')]
}


print('Brackets')
print(table(DT3$bracket))
print(table(DT3$bracketPC))

DT3[, use := unambiguousReads >= opt$cnt_thresh]
print('Use')
print(table(DT3$use))





DT1 <- dcast(DT1[name %in% c('primary', 'ambiguous'), ], sample ~ type, value.var = 'Count')
DT2 <- dcast(DT2[name %in% c('primary', 'ambiguous'), ], sample ~ type, value.var = 'PC')

DT <- merge(DT1, DT2, by='sample')
DT <- merge(DT, DT3[, .(sample, bracket, bracketPC, use)], by = 'sample')

fwrite(DT, file = file.path(opt$baseDir, opt$bbDir, 'bbsplit.stats.csv'), sep=',', col.names = TRUE)

pdf(file = file.path(opt$baseDir, opt$bbDir, 'bbsplit.stats.pdf'), width = 10, height = 10)
print(p1 / p2)
dev.off()
