#!/users/kimon.froussios/miniconda3/envs/mybasics/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
prefix <- args[2]             # for the collective PDF and HTML output files
ymax <- args[3]               # "NULL" for automatic
collectiveonly <- args[4]       # 'yes'/'no' useful to prevent file spam when many stats files are input
statsfiles <- args[4:length(args)]       # TSV input files: seq \t pos \t type \t count \t depth
# statsfiles <- '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/mutpe/process/B18_GCB/B18_GCB_KO_x25.extendedFrags_sorted_3-10.stats'

library(tidyverse)
library(ggrepel)
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

pdf(file.path(outdir, paste0(prefix, '.pdf')))

# sf <- statsfile[1]
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
  
  # Plot number of positions per mutation type
  print(posdata %>%
    filter((mutated)) %>%
    group_by(seq, type) %>%
    summarise(frequency = sum(count>0)) %>%
    ggplot(aes(x = type, y = frequency, fill=type)) + 
      facet_wrap( . ~ seq) +
      geom_bar(stat = 'identity', width = 0.9) +
      ylab('Number of positions') +
      ggtitle("Number of positions per mutation type.", subtitle=basename(sf))
  )
  
  # Plot number of reads per mutation type
  print(posdata %>%
    filter((mutated)) %>%
    group_by(seq, type) %>%
    summarise(frequency = sum(count)) %>%
    ggplot(aes(x = type, y = frequency, fill=type)) + 
      facet_wrap( . ~ seq) +
      geom_bar(stat = 'identity', width=0.9) +
      ylab('Number of reads') +
      ggtitle("Number of reads showing each mutation type.", subtitle=basename(sf))
  )
  
  # Pileup plot axis limits
  maxdepth <- max(posdata$depth)
  minpos <- min(posdata$pos)
  maxpos <- max(posdata$pos)
  maxfreq <- max(posdata$aggrfreq)
  if(ymax != 'NULL')
    maxfreq <- as.numeric(ymax)
  
  # Plot mutations pileup
  print(posdata %>% 
    filter(mutated) %>%
    select(seq, pos, aggrfreq, depth) %>%
    unique() %>%
    ggplot(aes(x = pos)) + 
      facet_wrap( . ~ seq) +
      geom_hline(aes(yintercept=quantile(aggrfreq, 0.25)), linetype='dotted', alpha=0.4) +
      geom_hline(aes(yintercept=median(aggrfreq)), linetype='dotdash', alpha=0.4) +
      geom_hline(aes(yintercept=quantile(aggrfreq, 0.75)), linetype='dashed', alpha=0.4) +
      geom_line(aes(y=depth / maxdepth * maxfreq * 1.1), colour='grey50') +
      geom_bar(aes(y = aggrfreq, fill=aggrfreq), stat = 'identity') +
      geom_label_repel(data=posdata %>% 
                       filter(mutated) %>%
                       select(pos, aggrfreq) %>%
                       unique() %>%
                       filter(aggrfreq >= quantile(posdata$aggrfreq, 0.95)),
                     aes(y=aggrfreq, label=pos), min.segment.length = 0, label.size = 0, fill="transparent", force=1) +
      scale_y_continuous('Frequency (bars)', sec.axis = sec_axis(trans= ~ . * maxdepth / maxfreq / 1.2, name = "Coverage (line)"), limits = c(0, 1)) +
      scale_fill_gradient(low='blue', high='magenta') +
      scale_x_continuous('Position', limits = c(minpos, maxpos)) +
      coord_cartesian(ylim=c(0, 1.1 * maxfreq)) +
      ggtitle("Fraction of reads with mutation at each position.", subtitle=basename(sf)) +
      theme(legend.position = 'none', 
            panel.background = element_rect(fill="white"), 
            panel.grid.major.y = element_line(colour='grey95'),
            panel.grid.minor.y = element_line(colour='grey95'))
  )
  
  if(collectiveonly=='yes'){
    # Plot interactive mutations pileup
    htmlwidgets::saveWidget(ggplotly(posdata %>% 
      filter(mutated) %>%
      mutate(type = factor(type, levels=c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G", "indel", "other"))) %>%
      ggplot(aes(x = pos, y = freq)) + 
        facet_wrap( . ~ seq) +
        geom_hline(aes(yintercept=quantile(aggrfreq, 0.25)), linetype='dotted', alpha=0.4) +
        geom_hline(aes(yintercept=median(aggrfreq)), linetype='dotdash', alpha=0.4) +
        geom_hline(aes(yintercept=quantile(aggrfreq, 0.75)), linetype='dashed', alpha=0.4) +
        geom_bar(aes(fill=type), stat = 'identity', position=position_stack()) +
        scale_x_continuous('Position', limits = c(minpos, maxpos)) +
        ylab("Frequency") +
        coord_cartesian(ylim=c(0, 1.1 * maxfreq)) +
        ggtitle("Fraction of reads with mutation at each position.", subtitle=basename(sf)) +
        guides(ncol=4, by.row=TRUE) +
        theme(panel.background = element_rect(fill="white"), 
              panel.grid.major.y = element_line(colour='grey95'),
              panel.grid.minor.y = element_line(colour='grey95'))
      ), 
      file.path(outdir, paste0(basename(sf), '_pileup_types.html'))
    )
  }
  
  # Concatenate all the files in one big table with the file as a variable, to plot overlay
  mega <- rbind(mega, posdata %>% filter(mutated) %>% mutate(file=basename(sf)))
}

dev.off()

# Plot mutations pileup
htmlwidgets::saveWidget(ggplotly(
  mega %>% 
    filter(mutated==TRUE) %>%
    select(seq, pos, aggrfreq, depth, file) %>%
    unique() %>%
    ggplot(aes(x = pos, y = aggrfreq, colour=file, fill=file)) + 
      facet_wrap( . ~ seq) +
      geom_bar(stat = 'identity', position=position_identity(), alpha=0.3) +
      scale_x_continuous("Position", limits = c(minpos, maxpos)) +
      scale_y_continuous('Coverage', limits = c(0, NA)) +
      ggtitle("Fraction of reads with mutation at each position.") +
      theme(panel.background = element_rect(fill="white"), 
            panel.grid.major.y = element_line(colour='grey95'),
        panel.grid.minor.y = element_line(colour='grey95'))
  ),
  file.path(outdir, paste0(prefix, '_pileups.html'))
)

# Plot mutations pileup
htmlwidgets::saveWidget(ggplotly(
  mega %>% 
    filter(mutated==TRUE) %>%
    select(seq, pos, aggrfreq, depth, file) %>%
    unique() %>%
    ggplot(aes(x = pos, y=depth, colour=file, fill=file)) + 
    facet_wrap( . ~ seq) +
    geom_line() +
    scale_x_continuous('Position', limits = c(minpos, maxpos)) +
    scale_y_continuous('Coverage', limits = c(0, NA)) +
    ggtitle("Fraction of reads with mutation at each position.") +
    theme(legend.position = 'bottom', 
          panel.background = element_rect(fill="white"), 
          panel.grid.major.y = element_line(colour='grey95'),
          panel.grid.minor.y = element_line(colour='grey95'))
  ),
  file.path(outdir, paste0(prefix, '_coverages.html'))
)

