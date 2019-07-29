#!/usr/bin/env Rscript

# Plots
# - mutation frequencies and read coverage statically for each input in a collective PDF
# - mutation frequencies interactively overlayed for all the input files
# - read coverage interactively overlayed for all the input files

## REQUIRES: pandoc
library(data.table)
library(dtplyr)
library(ggplot2)
library(ggrepel)
library(plotly)

#args <- c('~/', 'test', 'NULL', 'no', 'HDR1:280:6,HDR2:280:6', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process/in-vivo/HDR2/86754_B18_C2_G_e1_x25.extendedFrags_sorted_3-10.stats')
#args <- c('~/', 'test', 'NULL', 'no', '-:0:0', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process/in-vivo/HDRc/86752_B18_Cc_G_e1_x25.extendedFrags_sorted_1-2.stats')
#args <- c('~/', 'test', 'NULL', 'no', '-:0:0', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/yeap2015/process/pro/SRS1052796/SRR2229665_x25.extendedFrags_sorted_1-2.stats', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/yeap2015/process/pro/SRS1052796/SRR2229665_x25.extendedFrags_sorted_3-10.stats', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/yeap2015/process/pro/SRS1052796/SRR2229665_x25.extendedFrags_sorted11-30.stats')
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
prefix <- args[2]             # for the collective PDF and HTML output files
ymax <- args[3]               # "NULL" for automatic
collectiveonly <- args[4]     # 'yes'/'no' useful to prevent file spam when many stats files are input
trimmedstart <- 10            # shorten reference name. Keep from here
trimmedend <- 19              # shorten reference name. Keep to here
offsets <- args[5]            # Correct for deletions in references. "REF:START:LENGTH,REF:START:LENGTH" ie. "HDR2:280:6"
                              # For no corrections needed you can use "-:0:0". REF will be grepl'ed.
statsfiles <- args[6:length(args)]       # TSV input files: seq \t pos \t type \t count \t depth

# Expand the corrections string.
offsets <- strsplit(offsets, ',', fixed=TRUE)[[1]]
offsets <- as.data.frame(t(as.data.frame(strsplit(offsets, ':', fixed=TRUE))), stringsAsFactors=FALSE)
row.names(offsets) <- NULL

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

# sf <- statsfiles[1]
for (sf in statsfiles){
  # Input
  posdata <- fread(sf)
  if (dim(posdata)[1]==0) {
    warning(paste("No content found in ", sf))
    next
  }
  if (all(names(posdata) == c('amplicon', 'pos', 'type', 'freq', 'total')))
    names(posdata) <- c('seq', 'pos', 'type', 'count', 'depth')
  stopifnot(all(names(posdata) == c('seq', 'pos', 'type', 'count', 'depth')))

  # Correct for reference deletions
  for (ref in 1:dim(offsets)[1]) {
    #ref <- 1
    sel <- ( grepl(offsets[[ref,1]], posdata$seq) & 
             posdata$pos >= as.integer(offsets[[ref,2]]) )
    posdata$pos[sel] <- posdata$pos[sel] + as.integer(offsets[[ref,3]])
  }

  # Shorten the reference name so it fits in the facet bars.
  posdata[, seq := substr(seq, trimmedstart, trimmedend) ]

  # Create interruptions for coordinates without info.
  missing <- which(! 1:max(posdata$pos) %in% unique(posdata$pos))  # Maybe there is more info missing than just the deletions?
  if(length(missing) > 0){
    posdata <- rbind(posdata,
                     data.frame(seq=posdata$seq[1],      # Assumes a single reference per .stats file.
                                pos=missing,
                                type='ref_indel',
                                count=0,
                                depth=0) )
  }

  # Collapse and rename all insertion and deletions respectively.
  posdata[grepl('^[ACGT]$', type), type := 'ref']
  posdata[grepl('^+&', type), type := 'ins']
  posdata[grepl('^-$', type), type := 'del']
  
  # Fix order of categories.
  posdata[, type := ordered(type, levels=c("ref", "A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G", "ins", "del"))]
  
  # Calculate frequencies and aggregate frequencies by position
  posdata[, freq := count / depth]
  posdata[, mutated := freq != 0 & type != 'ref']
  posdata[(mutated), aggrfreq := sum(freq), by=c('seq', 'pos')]
  posdata[is.na(aggrfreq), aggrfreq := 0]

  # Plot number of positions per mutation type
  print(posdata %>%
    filter(mutated) %>%
    group_by(seq, type) %>%
    summarise(frequency = sum(count>0)) %>%
    ggplot(aes(x = type, y = frequency, fill=type)) +
      facet_grid( . ~ seq) +
      geom_bar(stat = 'identity', width = 0.9) +
      ylab('Number of positions') +
      ggtitle("Number of positions per mutation type.", subtitle=basename(sf)) +
      guides(fill="none") +
      theme(panel.background = element_rect(fill="white"), 
            panel.grid.major.y = element_line(colour="grey95"), 
            panel.grid.minor.y = element_line(colour="grey95"))
  )

  # Plot number of reads per mutation type
  print(posdata %>%
    filter(mutated) %>%
    group_by(seq, type) %>%
    summarise(frequency = sum(count)) %>%
    ggplot(aes(x = type, y = frequency, fill=type)) +
      facet_wrap( . ~ seq) +
      geom_bar(stat = 'identity', width=0.9) +
      #scale_y_log10() +
      ylab('Number of reads') +
      ggtitle("Number of reads showing each mutation type.", subtitle=basename(sf)) +
      guides(fill="none") +
      theme(panel.background = element_rect(fill="white"), 
            panel.grid.major.y = element_line(colour="grey95"), 
            panel.grid.minor.y = element_line(colour="grey95"))
  )

  # Pileup plot axis limits
  maxdepth <- max(posdata$depth)
  minpos <- min(posdata$pos)
  maxpos <- max(posdata$pos)
  maxfreq <- max(posdata %>% filter(type != 'ref') %>% select(aggrfreq))
  if(ymax != 'NULL')
    maxfreq <- as.numeric(ymax)
  # Automatically determine an additional zoomed level, to get around a few extremely talk peaks that squash everything else.
  zoomfreq <- quantile(posdata$aggrfreq, probs=c(0.99))
  
  
  # Plot mutations pileup
  pp <- posdata %>%
    filter(mutated) %>%
    select(seq, pos, aggrfreq, depth) %>%
    unique() %>%
    ggplot(aes(x = pos)) +
      facet_grid(seq ~ .) +
      geom_hline(aes(yintercept=quantile(aggrfreq, 0.25)), linetype='dotted', alpha=0.4) +
      geom_hline(aes(yintercept=median(aggrfreq)), linetype='dotdash', alpha=0.4) +
      geom_hline(aes(yintercept=quantile(aggrfreq, 0.75)), linetype='dashed', alpha=0.4) +
      geom_line(aes(y=depth / maxdepth * maxfreq * 1.1), colour='grey80') +
      geom_bar(aes(y = aggrfreq, fill=aggrfreq), stat = 'identity') +
      scale_y_continuous('Frequency (bars)', sec.axis = sec_axis(trans= ~ . * maxdepth / maxfreq / 1.2, name = "Coverage (line)")) +
      scale_fill_gradient(low='blue', high='magenta') +
      scale_x_continuous('Position', limits = c(minpos, maxpos), expand = c(0, 0)) +
      # coord_cartesian(ylim=c(NULL, 1.1 * maxfreq)) +
      ggtitle("Fraction of reads with mutation at each position.", subtitle=basename(sf)) +
      theme(legend.position = 'none',
            panel.background = element_rect(fill="white"),
            panel.grid.major.y = element_line(colour='grey95'),
            panel.grid.minor.y = element_line(colour='grey95'))
  print(pp + geom_label_repel(data=posdata %>%
                                filter(mutated) %>%
                                select(pos, aggrfreq) %>%
                                unique() %>%
                                filter(aggrfreq >= quantile(posdata$aggrfreq, 0.95)),
                              aes(y=aggrfreq, label=pos), 
                              force =1, min.segment.length = 0, label.size = 0, fill="transparent", xlim=c(minpos, maxpos))
  )
  print(pp + coord_cartesian(ylim=c(0, zoomfreq)) +
          geom_label_repel(data=posdata %>%
                             filter(mutated) %>%
                             select(pos, aggrfreq) %>%
                             unique() %>%
                             filter(aggrfreq >= quantile(posdata$aggrfreq, 0.95)),
                           aes(y=aggrfreq, label=pos), 
                           force = 1, min.segment.length = 0, label.size = 0, fill="transparent", xlim=c(minpos, maxpos), ylim=c(0, zoomfreq))
  )

  if(collectiveonly != 'yes'){
    # Plot interactive mutations pileup for each input
    htmlwidgets::saveWidget(ggplotly(posdata %>%
      filter(mutated) %>%
      mutate(type = factor(type, levels=c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G", "indel", "other"))) %>%
      ggplot(aes(x = pos, y = freq)) +
        facet_grid(seq ~ .) +
        geom_hline(aes(yintercept=quantile(aggrfreq, 0.25)), linetype='dotted', alpha=0.4) +
        geom_hline(aes(yintercept=median(aggrfreq)), linetype='dotdash', alpha=0.4) +
        geom_hline(aes(yintercept=quantile(aggrfreq, 0.75)), linetype='dashed', alpha=0.4) +
        geom_bar(aes(fill=type), stat = 'identity', position=position_stack()) +
        scale_x_continuous('Position', limits = c(minpos, maxpos)) +
        ylab("Frequency") +
        # coord_cartesian(ylim=c(NULL, 1.1 * maxfreq)) +
        ggtitle("Fraction of reads with mutation at each position.", subtitle=basename(sf)) +
        guides(ncol=4, by.row=TRUE) +
        theme(panel.background = element_rect(fill="white"),
              panel.grid.major.y = element_line(colour='grey95'),
              panel.grid.minor.y = element_line(colour='grey95'))
      ),
      file.path(outdir, paste0(basename(sf), '_pileup_MutTypes.html'))
    )
  }

  # Concatenate all the files in one big table with the file as a variable, to plot overlay
  mega <- rbind(mega, posdata %>% filter(mutated) %>% mutate(file=basename(sf)))
}

dev.off()

# Plot mutations collective pileups
htmlwidgets::saveWidget(ggplotly(
  mega %>%
    filter(mutated==TRUE) %>%
    select(seq, pos, aggrfreq, depth, file) %>%
    unique() %>%
    ggplot(aes(x = pos, y = aggrfreq, fill=file)) +
      facet_grid(seq ~ .) +
      geom_bar(stat = 'identity', position=position_identity(), alpha=0.4, colour='transparent') +
      scale_x_continuous("Position", limits = c(minpos, maxpos)) +
      scale_y_continuous('Frequency') +
      ggtitle("Fraction of reads with mutation at each position.") +
      theme(panel.background = element_rect(fill="white"),
            panel.grid.major.y = element_line(colour='grey95'),
        panel.grid.minor.y = element_line(colour='grey95'))
  ),
  file.path(outdir, paste0(prefix, '_pileups.html'))
)

# Plot collective coverages
htmlwidgets::saveWidget(ggplotly(
  mega %>%
    filter(mutated==TRUE) %>%
    select(seq, pos, aggrfreq, depth, file) %>%
    unique() %>%
    ggplot(aes(x = pos, y=depth, colour=file)) +
    facet_grid(seq ~ .) +
    geom_line(alpha=0.4) +
    scale_x_continuous('Position', limits = c(minpos, maxpos)) +
    scale_y_continuous('Coverage') +
    ggtitle("Fraction of reads with mutation at each position.") +
    theme(legend.position = 'bottom',
          panel.background = element_rect(fill="white"),
          panel.grid.major.y = element_line(colour='grey95'),
          panel.grid.minor.y = element_line(colour='grey95'))
  ),
  file.path(outdir, paste0(prefix, '_coverages.html'))
)
