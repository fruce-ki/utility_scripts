#!/users/kimon.froussios/miniconda3/envs/mybasics/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
prefix <- args[2]                       # for the collective PDF and HTML output files
ymax <- args[3]                         # "NULL" for automatic
statsfiles <- args[4:length(ags)]       # TSV input files: seq \t pos \t type \t count \t depth
# statsfiles <- '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/mutpe/process/B18_GCB/B18_GCB_KO_x0.extendedFrags_sorted_3-10.stats'

library(tidyverse)
library(plotly)

mega=data.frame(seq=NA_character_, 
                pos=NA_integer_, 
                type=NA_character_,
                count=NA_integer_,
                depth=NA_integer_,
                freq=NA_real_,
                mutated=NA,
                aggrfreq=NA_real_,
                file=NA_character_)

for (sf in statsfiles){
  # Input
  posdata <- read_tsv(sf) 
  if (! all(names(posdata) == c('amplicon', 'pos', 'type', 'freq', 'depth')))
    names(posdata) <- c('seq', 'pos', 'type', 'count', 'depth')
  stopifnot(all(names(posdata) == c('seq', 'pos', 'type', 'count', 'depth')))
  
  # Calculate frequencies and aggregate frequencies by position
  posdata <- posdata %>% 
    mutate(freq = count / depth, mutated = freq != 0) %>%
    group_by(seq, pos) %>%
    mutate(aggrfreq = sum(freq)) %>%
    ungroup()
  
  # Collapse all insertion and deletions respectively, and allow for unexpected categories
  deletions <- grepl('^del', posdata$type)
  insertions <- grepl('^ins', posdata$type)
  other <- grepl('^[^ACGT(ins)(del)(ref)]', posdata$type)
  posdata$type[deletions] <- 'indel' 
  posdata$type[insertions] <- 'indel'
  posdata$type[other] <- 'other'
