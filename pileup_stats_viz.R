#!~/miniconda3/envs/mybasics/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
collectivepdf <- args[2]
ymax <- args[3]                         # "NULL" for automatic
statsfiles <- args[4:length(ags)]       # seq pos type  count  depth
# statsfiles <- '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/mutpe/process/B18_GCB/B18_GCB_KO_x0.extendedFrags_sorted_3-10.stats'

library(tidyverse)
library(plotly)

pdf(file.path(outdir, collectivepdf))

# sf <- statsfile[1]
for (sf in statsfiles){
  # Input
  posdata <- read_tsv(sf) 
  # names(posdata) <- c('seq', 'pos', 'type', 'count', 'depth')
  fields <- names(posdata)
  stopifnot(all(fields == c('seq', 'pos', 'type', 'count', 'depth')))
  
  # Colculate frequencies and aggregate frequencies by position
  posdata <- posdata %>% 
    mutate(freq = count / depth, mutated = freq != 0) %>%
    group_by(seq, pos) %>%
    mutate(aggrfreq = sum(freq) / mean(depth)) %>%
    ungroup()
  
  # Collapse all insertion and deletions respectively, and allow for unexpected categories
  deletions <- grepl('^del', posdata$type)
  insertions <- grepl('^ins', posdata$type)
  other <- grepl('^[^ACGT(ins)(del)(ref)]', posdata$type)
  posdata$type[deletions] <- 'indel' 
  posdata$type[insertions] <- 'indel'
  posdata$type[other] <- 'other'
  
  mutPosByType <- posdata %>%
    filter((mutated)) %>%
    group_by(seq, type) %>%
    summarise(frequency = sum(count>0)) %>%
    ggplot(aes(x = type, y = frequency)) + 
      facet_wrap( . ~ seq) +
      geom_bar(stat = 'identity', width = 0.9) +
      ylab(paste0('Number of positions\n', basename(sf))) +
      ggtitle("Number of positions per mutation type.")
  
  mutReadsByType <- posdata %>%
    filter((mutated)) %>%
    group_by(seq, type) %>%
    summarise(frequency = sum(count)) %>%
    ggplot(aes(x = type, y = frequency)) + 
      facet_wrap( . ~ seq) +
      geom_bar(stat = 'identity', width=0.9) +
      ylab(paste0('Number of reads\n', basename(sf))) +
      ggtitle("Number of reads showing each mutation type.")
  
  maxdepth <- max(posdata$depth)
  minpos <- min(posdata$pos)
  maxpos <- max(posdata$pos)
  maxfreq <- max(posdata$freq)
  if(ymax != 'NULL')
    maxfreq <- as.numeric(ymax)

  mutTypeByPos <- posdata %>%
    filter((mutated)) %>%
    ggplot(aes(x = pos)) + 
      facet_wrap( . ~ seq) +
      geom_hline(aes(yintercept=quantile(freq, 0.25)), linetype='dotted', alpha=0.4) +
      geom_hline(aes(yintercept=median(freq)), linetype='dotdash', alpha=0.4) +
      geom_hline(aes(yintercept=quantile(freq, 0.75)), linetype='dashed', alpha=0.4) +
      geom_line(aes(y=depth / maxdepth * maxfreq * 1.1), colour='grey50') +
      geom_bar(aes(y = freq, fill = type), stat = 'identity', position = position_stack()) +
      scale_y_continuous(paste0('Frequency (bars)\n', basename(sf)), sec.axis = sec_axis(trans= ~ . * maxdepth / maxfreq / 1.2, name = "Coverage (line)"), limits = c(0, 1)) +
      scale_x_continuous('Position', limits = c(minx, maxx)) +
      ggtitle("Number of reads showing each mutation type.")
  
  mutTypeByPosZoom <- mutTypeByPos +
    coord_cartesian(ylim=c(0, 1.1 * maxfreq))
  
  print( mutPosByType )
  print( mutReadsByType )
  print( mutTypeByPosZoom )
  
  htmlwidgets::saveWidget(ggplotly(mutTypeByPosZoom), file.path(outdir, paste0(basename(sf), '.html')))
}

dev.off()





