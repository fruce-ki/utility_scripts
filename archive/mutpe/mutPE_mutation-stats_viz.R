#!/usr/bin/env Rscript

# Plots
# - mutation frequencies and read coverage statically for each input in a collective PDF
# - mutation frequencies interactively overlayed for all the input files
# - read coverage interactively overlayed for all the input files

## REQUIRES: pandoc
library(data.table)
library(ggplot2)
library(ggrepel)
library(plotly)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('~/', 'test', 'NULL', 'no', 'yes', 'HDR2:2301:6', 'B18', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process_vdj/in-vitro/HDR2/substractions/86740_B18_C2_s1_86732_B18_C2_r1.point.stats')
# args <- c('~/', 'test', 'NULL', 'no', 'no', 'HDR2:2301:6', 'B18', '/groups/pavri/Kimon/ursi/mutPEseq/round5/process_vdj/in-vitro/HDR2/substractions/86740_B18_C2_s1_86732_B18_C2_r1.point.stats')
# args <- c('~/', 'test', 'NULL', 'no', 'yes', '-:0:0', 'Ramos3', '/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/process_vdj/vdj3/56618_VDJ3_shRenilla.aln.point.stats')


trimmedstart <- 1             # shorten reference name. Keep from here
trimmedend <- 20              # shorten reference name. Keep to here
outdir <- args[1]             # output directory
prefix <- args[2]             # for the collective PDF and HTML output files
ymax <- args[3]               # "NULL" for automatic
collectiveonly <- args[4] == "yes"    # 'yes'/'no' useful to prevent file spam when many stats files are input
dobars <- args[5] == "yes"    # 'yes'/'no'. No will show only the pileup plots. Useful when individual mutation types are not tracked.
offsets <- args[6]            # Correct for deletions in references. "REF:START:LENGTH,REF:START:LENGTH" ie. "HDR2:280:6"
                              # For no corrections needed you can use "-:0:0". REF will be grepl'ed.
vdj <- args[7]                # select for predefined ranges
statsfiles <- args[8:length(args)]       # TSV input files: seq \t pos \t type \t count \t depth

covfilter <- 100    # min read count per position
cntfilter <- 30    # min mutated read count per peak

# View frames
if (vdj == 'B18' || vdj == 'B18_1') {
  showfrom <- 2022
  showto <- 2381
} else if (vdj == 'B18_1a') {
  showfrom <- 2022
  showto <- 2209
} else if (vdj == 'B18_1b') {
  showfrom <- 2195
  showto <- 2381
} else if (vdj == 'CH12' || vdj == 'CH12_1') {
  showfrom <- 1694
  showto <- 2059
} else if (vdj == 'CH12_1a') {
  showfrom <- 1694
  showto <- 1904
} else if (vdj == 'CH12_1b') {
  showfrom <- 1874
  showto <- 2077
} else if (vdj == 'Ramos1') {
  showfrom <- 3684
  showto <- 3929
} else if (vdj == 'Ramos2') {
  showfrom <- 3927
  showto <- 4147
} else if (vdj == 'Ramos3') {
  showfrom <- 4148
  showto <- 4375
} else if (vdj == 'Ramos') {
  showfrom <- 3684
  showto <- 4375
} else {
  stop(paste('Unknown view window', vdj, '.'))
}

# CDRs and hotspots
if (vdj == 'B18' || vdj == 'B18_1' || vdj == 'B18_1a' || vdj == 'B18_1b') {
  cdrs <- data.frame(xleft=c(2112, 2169, 2300) - 0.3,
                     xright=c(2126, 2216, 2353) + 0.3 )
  bed <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRCY_RGYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bed <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRCY_RGYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else if (vdj == 'CH12' || vdj == 'CH12_1' || vdj == 'CH12_1a' || vdj == 'CH12_1b') {
  cdrs <- data.frame(xleft=c(1784, 1841, 1988) - 0.3,
                     xright=c(1798, 1891, 2026) + 0.3 )
  bed <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRCY_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bed <- fread('/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRCY_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else if (grepl('Ramos', vdj)) {
  cdrs <- data.frame(xleft=c(3899, 3980, 4107) - 0.3,
                     xright=c(3943, 4021, 4181) + 0.3 )
  bed <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRCY_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bed <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRCY_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else {
  stop('No BEDfile defined for the viewing window', vdj, '.')
}
# Clean up BED
names(bed) <- c('chr', 'xleft', 'xright', 'seq')
bed[, xleft := xleft + 1]   # Convert coordinate ranges from 0-based right-open to 1-based right-closed
bed[, xleft := xleft - 0.3] # Pad segment ends to more closely line up with bar edges rather than centres
bed[, xrightt := xright + 0.3]
warmspots <- unique(bed[!grepl('[Aa][Gg][Cc][Tt]', seq, perl=TRUE), .(xleft, xright)])
hotspots <- unique(bed[grepl('[Aa][Gg][Cc][Tt]', seq, perl=TRUE), .(xleft, xright)]) # avoiding duplicate matches for palindromes

# Crashes are silent and mysterious because of the open PDF output file.
# So I have ot make checks before opening the PDF stream.
if (! dir.exists(outdir)) {
  stop(paste(outdir, "does not exist!"))
}
stopifnot(args[4] %in% c('yes', 'no'))
stopifnot(args[5] %in% c('yes', 'no'))
stopifnot(grepl('.*:\\d+:\\d+', offsets, perl=TRUE))
for (f in statsfiles) {
    if (! file.exists(f)) {
      stop(paste(f, "does not exist!"))
    }
}

# Expand the corrections string.
offsets <- strsplit(offsets, ',', fixed=TRUE)[[1]]
offsets <- as.data.frame(t(as.data.frame(strsplit(offsets, ':', fixed=TRUE))), stringsAsFactors=FALSE)
row.names(offsets) <- NULL

mega=data.table(chr=NA_character_,
                pos=NA_real_,
                type=NA_character_,
                count=NA_integer_,
                depth=NA_integer_,
                seq=NA_character_,
                newpos=NA_real_,
                mutated=NA,
                aggrcount=NA_real_,
                freq=NA_real_,
                aggrfreq=NA_real_,
                trusty=NA,
                aggrtrusty=NA,
                flag=NA,
                file=NA_character_)

# pdf('~/test.pdf', width=12, height=10)
pdf(file.path(outdir, paste0(prefix, '.pdf')), width=12, height=10)

counts <- data.table(file=statsfiles, max_coverage=NA_real_, min_coverage=NA_real_, mean_coverage=NA_real_, median_coverage=NA_real_)
setkey(counts, file)

# sf <- statsfiles[1]
for (sf in statsfiles){
  message(sf)
  # Is it an indel track? This affects how colour vectors must be defined.
  # The suffix pattern is given by the workflow script.
  isindel = grepl('.indel.stats', sf, fixed=TRUE)

  # Input
  info = file.info(sf)
  if (info$size == 0) {
    message(paste("No content found in ", sf))
    next
  }

  posdata <- fread(sf)
  if (dim(posdata)[1] == 0) {
    message(paste("No content found in ", sf))
    next
  }
  
  if (all(names(posdata) == c('amplicon', 'pos', 'type', 'freq', 'total')))
    names(posdata) <- c('chr', 'pos', 'type', 'count', 'depth')
  stopifnot(all(names(posdata) == c('chr', 'pos', 'type', 'count', 'depth')))
  
  # Correct coordinates for reference deletions
  posdata[, newpos := pos]
  for (ref in 1:dim(offsets)[1]) {
    #ref <- 1
    # Identify all rows from shift position onwards (on the relevant sequence)
    sel <- ( grepl(offsets[[ref,1]], posdata$chr) &
               posdata$pos >= as.integer(offsets[[ref,2]]) )
    # Create the shifted position. Retain the original position for the bedGraphs.
    posdata[(sel), newpos := pos + as.numeric(offsets[[ref,3]])]
  }
  
  # Apply the viewing window
  posdata <- posdata[newpos >= showfrom & newpos <= showto, ]
  
  # Don't even output graphs if there are no positions with OK coverage.
  if ( all(posdata$depth < cntfilter) ){
    message(paste("Insufficient coverage in", sf))
    next
  }
  
  # Shorten the reference name so it fits in the facet bars.
  posdata[, seq := regmatches(chr, regexpr('[A-Za-z0-9\\-]+_?[A-Za-z0-9\\-]*', chr, perl=TRUE))[1] ]
  
  # Collapse and rename all insertion and deletions respectively.
  posdata[grepl('^[NACGT]$', type), type := 'ref']
  if(!isindel){       # point mutations only or older style combined file
    posdata[grepl('+', type, fixed=TRUE), type := 'ins']
    posdata[grepl('-', type, fixed=TRUE), type := 'del']
  }

  
  # If there is nothing to show, don't try to. Not worth the code gymnastics to output blank plots.
  if(isindel){
    if(all(posdata[type == '-', count == 0])) {
      message(paste("No indels found in ", sf))
      next
    }
  } else {
    if(all(posdata$type == 'ref')) {
      message(paste("No mutations found in ", sf))
      next
    }
  }

  
  # Fix number, order and colour of categories for point mutations and combined files.
  levvec <- c("ref", "ins", "del", "any", "point_mut",
              "N>A", "C>A", "G>A", "T>A", 
              "N>C", "A>C", "G>C", "T>C", 
              "N>G", "A>G", "C>G", "T>G", 
              "N>T", "A>T", "C>T", "G>T")
  colvec <- c('ref'="transparent", 'ins'="grey50", 'del'="black",
              'any' = "aquamarine3", 'point_mut' = "cyan3",
              'N>A'="yellow", 'C>A'="gold", 'G>A'="orange", 'T>A'="darkorange",
              'N>C'="lightblue", 'A>C'="steelblue", 'G>C'="blue", 'T>C'="darkblue",
              'N>G'="pink", 'A>G'='violet', 'C>G'="magenta", 'T>G'="purple",
              'N>T'="indianred", 'A>T'="red", 'C>T'="red3", 'G>T'='red4')
  if(isindel){      # Don't know in advance which and how many categories there are.
    ipgen <- colorRampPalette(c('gold', 'darkred'))
    dpgen <- colorRampPalette(c('lightblue', 'darkblue'))
    v <- unique(posdata$type)
    i <- v[grepl('\\+[0-9]+', v, perl=TRUE)]
    i <- as.character(i[order(as.numeric(i))])
    d <- v[grepl('\\-[0-9]+', v, perl=TRUE)]
    d <- as.character(d[order(as.numeric(d), decreasing = TRUE)])
    levvec <- c('-', d, i)
    colvec <- c('grey75', dpgen(length(d)), ipgen(length(i)))
    names(colvec) <- levvec
  }
  posdata[, type := ordered(type, levels=levvec)]

  # Calculate frequencies and aggregate frequencies by position
  posdata[, mutated := count != 0 & (type != 'ref' | isindel)]       # There are no reference entries in the dedicated indel file
  posdata[(mutated), aggrcount := sum(count), by=c('seq', 'newpos')]    # Aggrecated mutation count (all mutation types per position)
  posdata[is.na(aggrcount), aggrcount := 0]
  posdata[(mutated), freq := count / depth]
  posdata[is.na(freq), freq := 0]
  posdata[, aggrfreq := aggrcount / depth]

  # Flip the frequencies of the deleted bases to negative Y, for better clarity
  if (isindel) {
    posdata[type=='-', freq := -freq]
    posdata[type=='-', aggrfreq := -aggrfreq]
  }
  
  # Trusty peaks
  posdata[, trusty := (mutated & abs(count) >= cntfilter & abs(depth) >= covfilter)]  # substraction tracks can have negative values
  posdata[grepl('N', type), trusty := FALSE]
  posdata[, aggrtrusty := aggrcount >= cntfilter, by=c("seq", "newpos")]  
  
  # Top 5% highest peaks
  posdata[(mutated), flag := aggrfreq >= quantile(aggrfreq, 0.95)]
  
  # Fill in gaps in the coordinates.
  setorder(posdata, seq, newpos)
  # Can't just look for missing positions between min and max, because there could be multiple sequences.
  # Identify the rows between which coordinates "jump". 
  missing <- data.table(from = which(c(vapply(1:(dim(posdata)[1]-1), function(i) {
                                  (posdata[i+1,newpos] > posdata[i, newpos+1]) & (posdata[i+1, seq] == posdata[i, seq])
                                }, logical(1)),
                                FALSE)), # Can't be anything missing after the last row.
                        to = which(c(FALSE,  # Can't be anything missing before the first row.
                               vapply(2:(dim(posdata)[1]), function(i) {
                                  (posdata[i-1,newpos] < posdata[i, newpos-1]) & (posdata[i-1, seq] == posdata[i, seq])
                                }, logical(1)) )) )
  if (nrow(missing) != 0) {
    for(i in 1:dim(missing)[1]){
      # How many rows to create
      n <- ceiling(posdata[missing[i, to], newpos] - posdata[missing[i, from], newpos] - 1)
      posdata <- rbind(posdata,
                       data.table(chr = rep(posdata[missing[i, from], chr], n),
                                  pos = NA_real_,
                                  type = NA_character_,
                                  count = 0,
                                  depth = 0,
                                  seq = rep(posdata[missing[i, from], seq], n),
                                  newpos = (posdata[missing[i, from], newpos] + 1):(posdata[missing[i, from], newpos] + n),
                                  mutated = FALSE,
                                  aggrcount = 0,
                                  freq = 0,
                                  aggrfreq = 0,
                                  trusty = FALSE,
                                  aggrtrusty = FALSE,
                                  flag = NA) )
    }
  }

  # Coverage statistics
  counts[sf, min_coverage := min(posdata$depth, na.rm=TRUE)]
  counts[sf, max_coverage := max(posdata$depth, na.rm=TRUE)]
  counts[sf, mean_coverage := mean(posdata$depth, na.rm=TRUE)]
  counts[sf, median_coverage := median(posdata$depth, na.rm=TRUE)]
  
  
  # Plot interactive mutations pileup. Stacked rate of each mutation type.
  if(!collectiveonly) {
    if (any(posdata$mutated)) {
      htmlwidgets::saveWidget(
        ggplotly(
          ggplot(posdata[(mutated), ],
                 aes(x = newpos, y = count, label=depth)) +
            facet_grid(seq ~ .) +
            geom_hline(aes(yintercept=quantile(aggrcount, 0.25)), linetype='dotted', alpha=0.4) +
            geom_hline(aes(yintercept=median(aggrcount)), linetype='dotdash', alpha=0.4) +
            geom_hline(aes(yintercept=quantile(aggrcount, 0.75)), linetype='dashed', alpha=0.4) +
            geom_bar(aes(fill=type), stat = 'identity', position=position_stack()) +
            scale_fill_manual(values=colvec) +
            # coord_cartesian(ylim=c(NULL, 1.1 * maxfreq)) +
            labs(title="Event rate per position", subtitle=basename(sf), y='Reads', x='Position') +
            guides(ncol=4, by.row=TRUE) +
            theme(panel.background = element_rect(fill="white"),
                  panel.grid.major.y = element_line(colour='grey95'),
                  panel.grid.minor.y = element_line(colour='grey95')),
          dynamicTicks = TRUE),
        file.path(outdir, paste0(basename(sf), '_pileup_MutTypes.html'))
      )
    }
  }
  
  # Plot number of positions per mutation type
  if (dobars){
    if (dim(posdata[mutated & type != '-',])[1] != 0) {
      print( # For an indel file, show only the indel declarations stats. The deleted bases are redundant and confusing here.
        
        ggplot(unique( posdata[mutated & type != '-', ] [, frequency := sum(count>0), by=.(seq, type)] [, .(seq, type, frequency)] ),
               aes(x = type, y = frequency, fill=type)) +
          facet_grid( . ~ seq) +
          geom_bar(stat = 'identity', width = 0.9) +
          scale_fill_manual(values=colvec) +
          ylab('Number of positions') +
          ggtitle("Positions per mutation type", subtitle=basename(sf)) +
          guides(fill="none") +
          theme(panel.background = element_rect(fill="white"),
                panel.grid.major.y = element_line(colour="grey95"),
                panel.grid.minor.y = element_line(colour="grey95"),
                axis.text.x = element_text(angle=60, hjust=1, vjust=1))
      )
  
      # Plot number of reads per mutation type
      print(
        ggplot(unique( posdata[mutated & type != '-', ] [, frequency := sum(count), by = .(seq, type)] [, .(seq, type, frequency)] ),
               aes(x = type, y = frequency, fill=type)) +
          facet_wrap( . ~ seq) +
          geom_bar(stat = 'identity', width=0.9) +
          scale_fill_manual(values=colvec) +
          # scale_y_sqrt() +
          ylab('Number of reads') +
          ggtitle("Reads per mutation type", subtitle=basename(sf)) +
          guides(fill="none") +
          theme(panel.background = element_rect(fill="white"),
                panel.grid.major.y = element_line(colour="grey95"),
                panel.grid.minor.y = element_line(colour="grey95"),
                axis.text.x = element_text(angle=60, hjust=1, vjust=1))
      )
    }
  }
  
  # next

  # Plot mutations pileup
  
  COV <- unique( posdata[, .(seq, newpos, depth)] )
  if (!isindel){
    DATA1 <- unique( posdata[(mutated & aggrfreq >= 0), .(seq, newpos, depth, aggrfreq, aggrtrusty)] )
    DATA2 <- unique( posdata[(mutated), .(seq, newpos, depth, aggrfreq, aggrtrusty)] )
  } else {
    DATA1 <- DATA2 <- unique( posdata[(mutated), .(seq, newpos, depth, aggrfreq, aggrtrusty)] )
  }
  minfreq <- min(DATA1$aggrfreq, na.rm=TRUE)
  maxfreq <- ifelse(any(DATA1$aggrfreq > 0),
                    max(DATA1[(aggrfreq > 0), aggrfreq], na.rm=TRUE),  # Use the largest positive value, if any
                    abs(max(DATA1[(aggrfreq >= 0), aggrfreq], na.rm=TRUE)) ) # use symmetric to the smallest negative value
  minpos <- min(posdata$newpos, na.rm=TRUE)
  maxpos <- max(posdata$newpos, na.rm=TRUE)
  maxdepth <- max(posdata$depth, na.rm=TRUE)
  # # Don't let low coverage positions force a zoom out
  # if (any(posdata$aggrtrusty & posdata$aggrfreq > 0, na.rm=TRUE)){
  #   maxfreq <- max(posdata[(aggrtrusty & aggrfreq > 0), aggrfreq], na.rm=TRUE)
  # } else if (any(posdata$mutated)) {
  #   # If there's nothing sufficiently reliable, use the st.dev. to exclude positions with outlier low coverage
  #   maxfreq = max(abs(posdata[depth >= (median(posdata$depth, na.rm=TRUE) - 3 * sd(posdata$depth, na.rm=TRUE)) & aggrfreq != 0, aggrfreq]), na.rm=TRUE)
  #   if(is.infinite(maxfreq)) {
  #     # If no mutations within 3SDs from median, just take whatever frequency there is.
  #     maxfreq = max(posdata[aggrfreq > 0, aggrfreq], na.rm=TRUE)
  #   }
  # } else {
  # Nothing to show anyway, so scale it to the smallest negative value (there have to be some values if we got this far).
  # maxfreq = abs(max(posdata$aggrfreq, na.rm=TRUE))
  # }
  # Unlike frequency, poor coverage can't create high count peaks.
  maxcount <- max(posdata[type != 'ref', abs(aggrcount)], na.rm=TRUE)
  # Manual zoom
  if(ymax != 'NULL') {
    maxfreq <- as.numeric(ymax)
  }
  # Automatically determine an additional zoomed level, to get around a few extremely tall peaks that squash everything else.
  #zoomfreq <- quantile(posdata$aggrfreq, probs=c(0.99))
  
  # Base
  pp <- ggplot(DATA2,
               aes(x = newpos, y = aggrfreq, fill = aggrtrusty)) +
    facet_grid(seq ~ .) +
    scale_x_continuous(expand=c(0,1)) +
    labs(y='Frequency', x="Position") +
    theme(legend.position = 'none',
          axis.text=element_text(size=9),
          axis.title=element_text(size=9),
          panel.background = element_rect(fill="white"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  
  # Three annotation styles.
  
  # zebra1 <- pp  # Stripe hotspots, dash CDRs below, zoomed in
  # zebra2 <- pp  # Stripe hotspots, dash CDRs above, zoomed out
  # dashy <- pp  # Stripe CDRs, dash hotspots below
  # 
  # if (nrow(warmspots) > 0) {
  #   for (i in 1:nrow(warmspots)) {
  #     zebra1 <- zebra1 +
  #       geom_rect(xmin=warmspots[i,]$xleft, xmax=warmspots[i,]$xright, ymin=0, ymax=1.5, fill='#FFEEBB')
  #   }
  #   for (i in 1:nrow(warmspots)) {
  #     zebra2 <- zebra2 +
  #       geom_rect(xmin=warmspots[i,]$xleft, xmax=warmspots[i,]$xright, ymin=-1, ymax=1.5, fill='#FFEEBB')
  #   }
  # }
  # if (nrow(hotspots) > 0) {
  #   for (i in 1:nrow(hotspots)) {
  #     zebra1 <- zebra1 +
  #       geom_rect(xmin=hotspots[i,]$xleft, xmax=hotspots[i,]$xright, ymin=0, ymax=1.5, fill='#FFD6BB')
  #   }
  #   for (i in 1:nrow(hotspots)) {
  #     zebra2 <- zebra2 +
  #       geom_rect(xmin=hotspots[i,]$xleft, xmax=hotspots[i,]$xright, ymin=-1, ymax=1.5, fill='#FFD6BB')
  #   }
  # }
  # if (nrow(cdrs) > 0) {
  #   zebra1 <- zebra1 +
  #     geom_segment(data=cdrs,
  #                  aes(x=xleft, xend=xright, y= -0.05 * maxfreq, yend= -0.05 * maxfreq), inherit.aes=FALSE,
  #                  colour='grey75', size=2)
  #   zebra2 <- zebra2 +
  #     geom_segment(data=cdrs,
  #                  aes(x=xleft, xend=xright, y= 1.1, yend= 1.1), inherit.aes=FALSE,
  #                  colour='grey50', size=2)
  #   for (i in 1:nrow(cdrs)) {
  #     dashy <- dashy +
  #       geom_rect(xmin=cdrs[i,]$xleft, xmax=cdrs[i,]$xright, ymin=-1, ymax=1.5, fill='grey95')
  #   }
  # }
  # # Keep dashes over the boxes, even if it means repeating code
  # if (nrow(warmspots) > 0) {
  #   dashy <- dashy +
  #     geom_segment(data=warmspots,
  #                  aes(x=xleft, xend=xright, y= -0.05 * maxfreq, yend= -0.05 * maxfreq), inherit.aes=FALSE,
  #                  colour='gold2', size=3)
  # }
  # if (nrow(hotspots) > 0) {
  #   dashy <- dashy +
  #     geom_segment(data=hotspots,
  #                aes(x=xleft, xend=xright, y= -0.05 * maxfreq, yend= -0.05 * maxfreq), inherit.aes=FALSE,
  #                colour='orangered', size=3)
  # }

  # Annotated version
  covscale1 <- 1.3  # Don't change the value for reuse with other plots. 
  pannot <- pp
  if (nrow(warmspots) > 0) {
    for (i in 1:nrow(warmspots)) {
      pannot <- pannot +
        geom_rect(xmin=warmspots[i,]$xleft, xmax=warmspots[i,]$xright, ymin=0, ymax=covscale1 * maxfreq, fill='#FFEEBB')
    }
  }
  if (nrow(hotspots) > 0) {
    for (i in 1:nrow(hotspots)) {
      pannot <- pannot +
        geom_rect(xmin=hotspots[i,]$xleft, xmax=hotspots[i,]$xright, ymin=0, ymax=covscale1 * maxfreq, fill='#FFD6BB')
    }
  }
  if (nrow(cdrs) > 0) {
    pannot <- pannot +
      geom_segment(data=cdrs,
                   aes(x=xleft, xend=xright, y= -0.05 * maxfreq, yend= -0.05 * maxfreq), inherit.aes=FALSE,
                   colour='grey50', size=2)
  }
  pannot <- pannot +
    geom_hline(yintercept=quantile(posdata[(mutated), aggrfreq], 0.25), linetype='dotted', alpha=0.3, color='forestgreen', size=0,5) +
    geom_hline(yintercept=median(posdata[(mutated), aggrfreq]), linetype='dotdash', alpha=0.3, color='forestgreen', size=0,5) +
    geom_hline(yintercept=quantile(posdata[(mutated), aggrfreq], 0.75), linetype='dashed', alpha=0.3, color='forestgreen', size=0,5) +
    geom_line(data=COV,
              aes(x=newpos, y=depth / maxdepth * maxfreq * covscale1), inherit.aes=FALSE, colour='steelblue1', size=0.25) +
    geom_bar(data=DATA1, stat = 'identity', width=0.9) +
    geom_hline(yintercept=0, size=0.25) +
    geom_label_repel(data=unique(posdata[(flag), .(newpos, aggrfreq)]), 
                     aes(x=newpos, y=aggrfreq, label=newpos), inherit.aes=FALSE,
                     min.segment.length = 0, label.size = 0, fill="transparent", colour='red',
                     xlim=c(minpos, maxpos), ylim=c(minfreq, maxfreq * covscale1 * 0.98)) +
    scale_y_continuous(sec.axis = sec_axis(trans= ~ . * maxdepth / maxfreq / covscale1, name = "Coverage"), expand=c(0,0)) +
    coord_cartesian(xlim=c(showfrom, showto),
                    ylim=c(ifelse(isindel, minfreq, -0.1 * maxfreq), maxfreq * covscale1 * 1.02)) +
    labs(title="Event rate per position", subtitle=basename(sf)) +
    scale_fill_manual(values=c('FALSE'="royalblue", 'TRUE'="black"), na.value="transparent")
  
  # Simpler versions
  covscale2 <- 1.1
  pplain <- pp
  if (nrow(cdrs) > 0) {
    for (i in 1:nrow(cdrs)) {
      pplain <- pplain +
        geom_rect(xmin=cdrs[i,]$xleft, xmax=cdrs[i,]$xright, ymin=ifelse(isindel, minfreq, -0.1 * maxfreq), ymax=covscale2 * maxfreq, fill='grey95')
    }
  }
  if (nrow(warmspots) > 0) {
    pplain <- pplain +
      geom_segment(data=warmspots,
                   aes(x=xleft, xend=xright, y= ifelse(isindel, minfreq, -0.05 * maxfreq), yend= ifelse(isindel, minfreq, -0.05 * maxfreq)), inherit.aes=FALSE,
                   colour='gold2', size=3)
  }
  if (nrow(hotspots) > 0) {
    pplain <- pplain +
      geom_segment(data=hotspots,
                   aes(x=xleft, xend=xright, y= ifelse(isindel, minfreq, -0.05 * maxfreq), yend= ifelse(isindel, minfreq, -0.05 * maxfreq)), inherit.aes=FALSE,
                   colour='orangered', size=3)
  }
  
  pplain1 <- pplain +
    geom_line(data=COV,
              aes(x=newpos, y=depth / maxdepth * maxfreq * covscale2), inherit.aes=FALSE, colour='steelblue1', size=0.25) +
    geom_bar(data=DATA1, stat = 'identity', width=0.9) +
    geom_hline(yintercept=0, size=0.6) +
    scale_y_continuous(sec.axis = sec_axis(trans= ~ . * maxdepth / maxfreq / covscale2, name = "Coverage"), expand=c(0,0)) +
    coord_cartesian(xlim=c(showfrom, showto),
                    ylim=c(ifelse(isindel, minfreq, -0.1 * maxfreq), maxfreq * covscale2 * 1.02)) +
    scale_fill_manual(values=c('FALSE'="royalblue", 'TRUE'="black"), na.value="transparent")
  pplain2 <- pplain +
    geom_line(data=COV,
              aes(x=newpos, y=depth / maxdepth * maxfreq * covscale2), inherit.aes=FALSE, colour='steelblue1', size=0.25) +
    geom_bar(data=DATA1, stat = 'identity', width=0.9) +
    geom_hline(yintercept=0, size=0.6) +
    scale_y_continuous(sec.axis = sec_axis(trans= ~ . * maxdepth / maxfreq / covscale2, name = "Coverage"), expand=c(0,0)) +
    coord_cartesian(xlim=c(showfrom, showto),
                    ylim=c(ifelse(isindel, minfreq, -0.1 * maxfreq), maxfreq * covscale2)) +
    scale_fill_manual(values=c('FALSE'="black", 'TRUE'="black"), na.value="transparent")
  pplain3 <- pplain +
    geom_bar(data=DATA1, stat = 'identity', width=0.9) +
    geom_hline(yintercept=0, size=0.6) +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(xlim=c(showfrom, showto),
                    ylim=c(ifelse(isindel, minfreq, -0.1 * maxfreq), maxfreq * covscale2)) +
    scale_fill_manual(values=c('FALSE'="black", 'TRUE'="black"), na.value="transparent")
  
  # Zoomed out version
  pzoom <- pp
  if (nrow(warmspots) > 0) {
    for (i in 1:nrow(warmspots)) {
      pzoom <- pzoom +
        geom_rect(xmin=warmspots[i,]$xleft, xmax=warmspots[i,]$xright, ymin=-1, ymax=1.1, fill='#FFEEBB')
    }
  }
  if (nrow(hotspots) > 0) {
    for (i in 1:nrow(hotspots)) {
      pzoom <- pzoom +
        geom_rect(xmin=hotspots[i,]$xleft, xmax=hotspots[i,]$xright, ymin=-1, ymax=1.1, fill='#FFD6BB')
    }
  }
  if (nrow(cdrs) > 0) {
    pzoom <- pzoom +
      geom_segment(data=cdrs,
                   aes(x=xleft, xend=xright, y= 1.1, yend= 1.1), inherit.aes=FALSE,
                   colour='grey50', size=2)
  }
  pzoom <- pzoom + 
    geom_line(data=COV,
              aes(x=newpos, y=depth / maxdepth), inherit.aes=FALSE, colour='steelblue1', size=0.25) +
    geom_bar(data=DATA2, stat = 'identity', width=0.9) +
    geom_hline(yintercept=0, size=0.25) +
    scale_y_continuous(sec.axis = sec_axis(trans= ~ . * maxdepth, name = "Coverage"), expand=c(0,0)) +
    scale_fill_manual(values=c('FALSE'="royalblue", 'TRUE'="black"), na.value="transparent") +
    coord_cartesian(xlim=c(showfrom, showto), ylim=c(-1, 1.15))
  
  print( pannot / pplain1 / pplain2 / pplain3 / pzoom )
  
  # De-factorise apples and oranges so I can bind point mutation and indel tables.
  posdata <- posdata[, type := as.character(type)]

  
  # Create bedGraph of mutations
  # The decimal offset coordinates for indel declarations might freak out a genome browser and I don't need to solve this (yet).
  # Use the normal mapped coordinates for the bedGraph, not the shifted ones used for peak labelling in the plots.
  if (!isindel){  
    # bG <- rbind(unique(posdata[(mutated), .(chr, pos-1, pos, aggrcount * 1)] [(V4 < 0), V4 := 0]),  
    #           unique(posdata[(mutated), .(chr, pos-1, pos, depth/maxdepth * maxcount * -1)]) )  # scale to counts and flip to negatives
    # multiplying by 1 to force column rename and make it compatible with the second table.
    # if the track is a substraction track, negative frequencies will clash with the coverage info I want to put as negative values, so remove them.
    # OR don't overcomplicate the track:
    bG <- unique(posdata[(mutated), .(chr, pos-1, pos, aggrcount)] [(aggrcount < 0), aggrcount := 0])
    if (dim(bG)[1] != 0) {
      fwrite(bG, file=sub('stats', paste0('cov', maxdepth, '.bedGraph'), sf), col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)
    } else {
      message("No mutations to put in the bedGraph.")
    }
  }
  
  # Concatenate all the files in one big table with the file as a variable, to plot overlay
  mega <- rbind(mega,
                posdata[(mutated), ] [, file := basename(sf)] )
}

dev.off()


# Interactive pileup overlay and coverage overlay.
mega <- mega[complete.cases(mega),]   # mainly to get rid of the dummy row
if (! dim(mega[(mutated), ])[1] == 0) {
  # Pileup plot axis limits
  maxdepth <- max(mega$depth)
  minpos <- min(mega$pos)
  maxpos <- max(mega$pos)
  maxfreq <- max(mega[type != 'ref', aggrfreq], na.rm=TRUE)
  
  if(ymax != 'NULL')
    maxfreq <- as.numeric(ymax)
  
  # Plot mutations collective pileups
  if(any(mega$aggrtrusty)){
    htmlwidgets::saveWidget(
      ggplotly(
        ggplot(unique( mega[(aggrtrusty), .(seq, pos, aggrfreq, depth, file)] ), 
               aes(x = pos, y = aggrfreq, fill=file, label=depth)) +
          facet_grid(seq ~ .) +
          geom_bar(stat = 'identity', position=position_identity(), alpha=0.4, colour='transparent') +
          # scale_x_continuous(limits = c(minpos, maxpos)) +
          labs(title=paste0("Events with Nreads>=", cntfilter," and depth>=", covfilter), y='Frequency', x='Position') +
          theme(panel.background = element_rect(fill="white"),
                panel.grid.major.y = element_line(colour='grey95'),
                panel.grid.minor.y = element_line(colour='grey95')),
        dynamicTicks = TRUE),
      file.path(outdir, paste0(prefix, '_pileups.html'))
    )
  } else {
    message("Nothing reliable to output to the collective pileup.")
  }
  
  # Plot collective coverages
  htmlwidgets::saveWidget(
    ggplotly(
      ggplot(mega, aes(x = pos, y=depth, colour=file)) +
        facet_grid(seq ~ .) +
        geom_line(alpha=0.4) +
        coord_cartesian(ylim=c(0, maxdepth)) +
        # scale_x_continuous('Position', limits = c(minpos, maxpos)) +
        labs(title = "Position coverage", y="Reads", x= "Position") +
        theme(panel.background = element_rect(fill="white"),
              panel.grid.major.y = element_line(colour='grey95'),
              panel.grid.minor.y = element_line(colour='grey95')),
      dynamicTicks = TRUE),
    file.path(outdir, paste0(prefix, '_coverages.html'))
  )
}

fwrite(counts, file=file.path(outdir, paste0(prefix, '_coverages.tsv')), sep="\t")

# warnings()
