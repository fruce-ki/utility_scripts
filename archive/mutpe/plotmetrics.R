#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)

options(scipen = 100)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/process_vdj/vdj1/metrics.txt', "~/test.pdf", "1-1", "2-2", "3-3", "4-4", "5-15", "16-30")
metrics <- args[1]
out <- args[2]
strata <- args[3:length(args)]

flow <- c("fastq R1", "fastq R2", "trim3 R1", "trim3 R2", "trim53 R1", "trim53 R2", "extendedFrags", "extendedFrags_fltr", "aln",
          paste0('aln_', strata))
flow2 <- c("Input R1", "Input R2", "Quality trim R1", "Quality trim R2", "Adaptor trim R1", "Adaptor trim R2", "Merge", "Length", "Align", strata)

DT <- unique(fread(metrics))
stopifnot(dim(DT)[1] > 0)

# Analyze the file names to create grouping factors
DT[, file := sub('extendedFrags\\.fltr', 'extendedFrags_fltr', file)]
DT[, file := sub('aln\\.(?!bam)', 'aln_', file, perl=TRUE)]
DT[, c("sample", "stage") := tstrsplit(file, '.', fixed=TRUE, keep=c(1,2))]
DT[grepl('_1\\.fastq|_1\\.trim', file, perl=TRUE), stage := paste(stage, 'R1')]
DT[grepl('_2\\.fastq|_2\\.trim', file, perl=TRUE), stage := paste(stage, 'R2')]
DT[grepl('_[12]\\.fastq|_[12]\\.trim', file, perl=TRUE), sample := sub('_[12]$', '', sample)]
DT[, stage := factor(stage, ordered=TRUE, levels=flow)]
setorder(DT, sample, stage)

stopifnot(all(complete.cases(DT)))  # Probably changed strata definitions in workflow.sh and didn't update the call to this script.

# Calculate percentages
DT[, ref := max(count), by=sample]
DT[, fraction := count / ref]


p1 <- ggplot(DT, aes(x=stage, y=count, colour=sample, group=sample, label=substr(sample, 1, 30))) +
  geom_line(alpha=0.4, size=0.6) +
  geom_text_repel(data=DT[stage=="aln",], direction = "y", nudge_x = 5, segment.colour = 'black', segment.size = 0.3, segment.alpha = 0.6) +
  scale_x_discrete(labels=flow2) +
  theme_bw() +
  labs(title=sub("^.*mutPEseq/", "", dirname(metrics)), x="", y="Read count") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        legend.position="none")

p2 <- ggplot(DT, aes(x=stage, y=count, colour=sample, group=sample, label=substr(sample, 1, 30))) +
  geom_line(alpha=0.4, size=0.6) +
  geom_text_repel(data=DT[stage=="aln",], direction = "y", nudge_x = -7, nudge_y = -2, segment.colour = 'black', segment.size = 0.3, segment.alpha = 0.6) +
  scale_x_discrete(labels=flow2) +
  scale_y_log10() +
  annotation_logticks(sides = "rl") +
  labs(x="", y="") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        panel.grid.minor.y = element_blank(),
        legend.position="none")


p3 <- ggplot(DT, aes(x=stage, y=fraction, colour=sample, group=sample, label=substr(sample, 1, 30))) +
  geom_line(alpha=0.4, size=0.6) +
  geom_text_repel(data=DT[stage=="aln",], direction = "y", nudge_x = 5, segment.colour = 'black', segment.size = 0.3, segment.alpha = 0.6) +
  scale_x_discrete(labels=flow2) +
  labs(y="Count fraction") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        legend.position="none")

p4 <- ggplot(DT, aes(x=stage, y=fraction, colour=sample, group=sample, label=substr(sample, 1, 30))) +
  geom_line(alpha=0.4, size=0.6) +
  geom_text_repel(data=DT[stage=="aln",], direction = "y", nudge_x = -7, nudge_y = -2, segment.colour = 'black', segment.size = 0.3, segment.alpha = 0.6) +
  scale_x_discrete(labels=flow2) +
  scale_y_log10() +
  annotation_logticks(sides = "rl") +
  labs(y="") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        panel.grid.minor.y = element_blank(),
        legend.position="none")

# Save plot
pdf(file=out, paper="a4r", width=11, height=8 )
p1 + p2 + p3 + p4 
dev.off()

# Save table in people-friendly wide format
DT <- dcast(DT, sample ~ stage, value.var=c("count", "fraction"))
names(DT) <- sub(" ", "_", names(DT))
fwrite(DT, 
       file=sub(".pdf", ".tsv", out), 
       sep="\t")
