#!/usr/bin/env Rscript


## REQUIRES: pandoc
library(data.table)
library(ggplot2)
library(ggrepel)
library(plotly)
library(patchwork)


args <- commandArgs(trailingOnly = TRUE)
# args <- c('~', 'test', '-:0:0', 'Ramos', '0.2', '/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/process_vdj/mergers/56606-VDJ1_56612-VDJ2_56618-VDJ3_shRenilla.aln.16-30.point.stats')
# args <- c('~', 'test', '-:0:0', 'B18', '1', '/Volumes/groups/pavri/Kimon/ursi/mutPEseq/round5/process_paper/substractions/86745_B18_AIDER_s1_86736_B18_AIDKO_s1.point.stats')


trimmedstart <- 1             # shorten reference name. Keep from here
trimmedend <- 20              # shorten reference name. Keep to here
outdir <- args[1]             # output directory
prefix <- args[2]             # for the collective PDF and HTML output files
offsets <- args[3]            # Correct for deletions in references. "REF:START:LENGTH,REF:START:LENGTH" ie. "HDR2:280:6"
                              # For no corrections needed you can use "-:0:0". REF will be grepl'ed.
vdj <- args[4]                # select for predefined ranges
allelecutoff <- as.numeric(args[5])        # Ignore mutations greater than 30% of coverage.
statsfiles <- args[6:length(args)]       # TSV input files: seq \t pos \t type \t count \t depth

# Crashes are silent and mysterious because of the open PDF output file.
# So I have ot make checks before opening the PDF stream.
if (! dir.exists(outdir)) {
  stop(paste(outdir, "does not exist!"))
}
stopifnot(grepl('.*:\\d+:\\d+', offsets, perl=TRUE))
for (f in statsfiles) {
  if (! file.exists(f)) {
    stop(paste(f, "does not exist!"))
  }
}



# View range
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
  bedfw <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRCY.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrv <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_RGYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedfws <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRC.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrvs <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_GYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedfwn <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_SYC.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrvn <- fread('/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_GRS.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfw <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRCY.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrv <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_RGYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfws <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_WRC.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrvs <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_GYW.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfwn <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_SYC.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrvn <- fread('/Volumes/groups/pavri/Kimon/ursi/PRO/aux/VDJ_locus/B18_all_GRS.bed', skip=1)[V1=='B1-8hi_mm10_190516_181014-most-recent' & V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else if (vdj == 'CH12' || vdj == 'CH12_1' || vdj == 'CH12_1a' || vdj == 'CH12_1b') {
  cdrs <- data.frame(xleft=c(1784, 1841, 1988) - 0.3,
                     xright=c(1798, 1891, 2026) + 0.3 )
  bedfw <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRCY.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrv <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedfws <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrvs <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_GYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedfwn <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrvn <- fread('/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfw <- fread('/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRCY.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrv <- fread('/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfws <- fread('/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_WRC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrvs <- fread('/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_GYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfwn <- fread('/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrvn <- fread('/Volumes/groups/pavri/Kimon/ursi/CH12_MutPE/aux/VDJ_locus/CH12_VDJ_200122_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else if (grepl('Ramos', vdj)) {
  cdrs <- data.frame(xleft=c(3899, 3980, 4107) - 0.3,
                     xright=c(3943, 4021, 4181) + 0.3 )
  bedfw <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRCY.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrv <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedfws <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrvs <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_GYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedfwn <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  bedrvn <- fread('/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfw <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRCY.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrv <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_RGYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfws <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_WRC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrvs <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_GYW.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedfwn <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_SYC.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
  # bedrvn <- fread('/Volumes/groups/pavri/Kimon/ursi/Ramos_MutPE/aux/VDJ/Ramos_IgH_GRS.bed', skip=1)[V3 >= showfrom & V2 <= showto, .(V1, V2, V3, V4)]
} else {
  stop('No BEDfile defined for the viewing window', vdj, '.')
}
# Clean up BED
names(bedfw) <- c('chr', 'xleft', 'xright', 'seq')
names(bedrv) <- c('chr', 'xleft', 'xright', 'seq')
names(bedfws) <- c('chr', 'xleft', 'xright', 'seq')
names(bedrvs) <- c('chr', 'xleft', 'xright', 'seq')
names(bedfwn) <- c('chr', 'xleft', 'xright', 'seq')
names(bedrvn) <- c('chr', 'xleft', 'xright', 'seq')
# Isolate core coordinate of the pattern
bedfw[, x := xleft + 3]    # wrCy in 1-based coordinate from 0 bases coordinate range
bedrv[, x := xleft + 2]    # rGyw in 1-based coordinate
bedfws[, x := xleft + 3]   # wrC in 1-based coordinate
bedrvs[, x := xleft + 1]   # Gyw in 1-based coordinate
bedfwn[, x := NA_real_]
bedrvn[, x := NA_real_]

# Combine
bedfw[, fw := TRUE]
bedrv[, fw := FALSE]
bedfws[, fw := TRUE]
bedrvs[, fw := FALSE]
bedfwn[, fw := TRUE]
bedrvn[, fw := FALSE]
bedfw[, origin := 'long']
bedrv[, origin := 'long']
bedfws[, origin := 'short']
bedrvs[, origin := 'short']
bedfwn[, origin := 'cold']
bedrvn[, origin := 'cold']
bed <- rbind(bedfw, bedrv, bedfws[!x %in% bedfw$x,], bedrvs[!x %in% bedrv$x,], bedfwn, bedrvn) # remove overlaps (same hot C/G base)
setorder(bed, xleft)

# Convert coordinate ranges from 0-based right-open to 1-based right-closed
bed[, xleft := xleft + 1] 
# Pad segment ends to more closely line up with bar edges rather than centres
bed[, xleft := xleft - 0.3]
bed[, xright := xright + 0.3]
# Distinguish between strict and loose pattern. And lose overlaps between the short and long patterns.
# The seq field is already rev'comp'ed for reverse strand matches, by the way the BEDfiles were made.
hotspotslong <- bed[origin == 'long' & grepl('[Aa][Gg][Cc][Tt]', seq, perl=TRUE), ] [, origin := 'hotlong']
hotspotshort <- bed[origin == 'short' & grepl('[Aa][Gg][Cc]', seq, perl=TRUE) & !(xleft %in% hotspotslong$xleft | xright %in% hotspotslong$xright), ] [, origin := 'hotshort']
warmspotslong <- bed[origin == 'long' & !grepl('[Aa][Gg][Cc][Tt]', seq, perl=TRUE), ] [, origin := 'warmlong']
warmspotshort <- bed[origin == 'short' & !grepl('[Aa][Gg][Cc]', seq, perl=TRUE) & !(xleft %in% warmspotslong$xleft | xright %in% warmspotslong$xright), ] [, origin := 'warmshort']
coldspots <- bed[origin=='cold', ]
# Update designations.
bed <- rbind(hotspotslong, hotspotshort, warmspotslong, warmspotshort, coldspots)
# Prepare annotations. Lose duplicates from palindromes.
hotspotslong <- unique(bed[origin == 'hotlong', .(xleft, xright)])
hotspotshort <- unique(bed[origin == 'hotshort', .(xleft, xright)])
warmspotslong <- unique(bed[origin == 'warmlong', .(xleft, xright)])
warmspotshort <- unique(bed[origin == 'warmshort', .(xleft, xright)])
coldspots <- unique(bed[origin=='cold', .(xleft, xright)])

# Shift annotation coordinates so x axis starts at 1.
cdrs <- cdrs - showfrom + 1
warmspotslong <- warmspotslong - showfrom + 1
hotspotslong <- hotspotslong - showfrom + 1
warmspotshort <- warmspotshort - showfrom + 1
hotspotshort <- hotspotshort - showfrom + 1
coldspots <- coldspots - showfrom + 1


# Expand the coordinate corrections string.
offsets <- strsplit(offsets, ',', fixed=TRUE)[[1]]
offsets <- as.data.frame(t(as.data.frame(strsplit(offsets, ':', fixed=TRUE))), stringsAsFactors=FALSE)
row.names(offsets) <- NULL

# Outputs
# pdf('~/test.pdf', width=12, height=10)
pdf(file.path(outdir, paste0(prefix, '.pdf')), width=12, height=3)
mtxf <- file.path(outdir, paste0(prefix, '_matrices.txt'))
cat('', file=mtxf)


for (sf in statsfiles){
  # sf <- statsfiles[1]
  message(sf)
  info = file.info(sf)
  if (info$size == 0) {
    message(paste("No content found in ", sf))
    next
  }

  # Input
  posdata <- fread(sf)
  if (dim(posdata)[1] == 0) {
    message(paste("No content found in ", sf))
    next
  }
  if (all(names(posdata) == c('amplicon', 'pos', 'type', 'freq', 'total')))
    names(posdata) <- c('chr', 'pos', 'type', 'count', 'depth')
  stopifnot(all(names(posdata) == c('chr', 'pos', 'type', 'count', 'depth')))
  
  # Correct coordinates for reference deletions
  posdata[, newpos := as.numeric(pos)]
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
  
  # Shorten the reference name so it fits in the facet bars.
  posdata[, seq := regmatches(chr, regexpr('[A-Za-z0-9\\-]+_?[A-Za-z0-9\\-]*', chr, perl=TRUE))[1] ]
  
  # Collapse and rename all insertion and deletions respectively.
  posdata[grepl('^[NACGT]$', type), type := 'ref']
  posdata[grepl('+', type, fixed=TRUE), type := 'ins']
  posdata[grepl('-', type, fixed=TRUE), type := 'del']
  
  # If there is nothing to show, don't try to. Not worth the code gymnastics to output blank plots.
  if(all(posdata$type == 'ref')) {
    message(paste("No mutations found in ", sf))
    next
  }
  
  # Categorize mutations.
  posdata[, mutated := grepl('>', type, fixed=TRUE)]
  # by source base type
  posdata[(mutated), grouptype := 'A:T']
  posdata[mutated & (grepl('C>', type, fixed=TRUE) | grepl('G>', type, fixed=TRUE)), grouptype := 'G:C']
  # by presence on a warm or hot spot
  posdata[(mutated), sweetspotlong := FALSE]
  posdata[mutated & grepl('C>', type, fixed=TRUE) & pos %in% bed[fw & origin %in% c('hotlong', 'warmlong'), x], sweetspotlong := TRUE]
  posdata[mutated & grepl('G>', type, fixed=TRUE) & pos %in% bed[(!fw) & origin %in% c('hotlong', 'warmlong'), x], sweetspotlong := TRUE]
  posdata[(mutated), sweetspotshort := FALSE]
  posdata[mutated & grepl('C>', type, fixed=TRUE) & pos %in% bed[fw & origin %in% c('hotshort', 'warmshort'), x], sweetspotshort := TRUE]
  posdata[mutated & grepl('G>', type, fixed=TRUE) & pos %in% bed[(!fw) & origin %in% c('hotshort', 'warmshort'), x], sweetspotshort := TRUE]
  posdata[(mutated), sourspot := FALSE]
  posdata[mutated & grepl('C>', type, fixed=TRUE), 
          sourspot := vapply(posdata[mutated & grepl('C>', type, fixed=TRUE), pos], 
                             function(p) {
                               any(bed[fw & origin == 'cold', p > xleft & p < xright]) 
                             }, logical(1)) ] # Remember, the x have been padded with a decimal value
  posdata[mutated & grepl('G>', type, fixed=TRUE), 
          sourspot := vapply(posdata[mutated & grepl('G>', type, fixed=TRUE), pos], 
                             function(p) {
                               any(bed[(!fw) & origin == 'cold', p > xleft & p < xright]) 
                             }, logical(1)) ] # Remember, the x have been padded with a decimal value
  posdata[, mutcat := grouptype]
  posdata[grouptype == 'G:C', mutcat := 'other C:G']
  posdata[(sweetspotlong), mutcat := 'C:G in WRCY/RGYW']
  posdata[(sweetspotshort), mutcat := 'C:G in other WRC/GYW']
  posdata[(sourspot) & !(sweetspotlong | sweetspotshort), mutcat := 'C:G in SYC/GRS']      # not claimed by a hot/warm spot
  posdata[(sourspot) & (sweetspotlong | sweetspotshort), mutcat := 'C:G in WRC/GYW and SYC/GRS']       # conflict with hot/warm spot

  # Fix number, order and colour of categories for point mutations and combined files.
  levvec <- c("ref", "ins", "del", "any", "point_mut",
              "N>A", "N>T", "N>G", "N>C", 
                     "A>T", "A>G", "A>C", "A>N", 
              "T>A",        "T>G", "T>C", "T>N",
              "G>A", "G>T",        "G>C", "G>N", 
              "C>A", "C>T", "C>G",        "C>N")
  colvec <- c("ref"='black', "ins"='grey20', "del"='grey10', "any"='grey50', "point_mut"='grey40',
              "N>A"='grey90', "N>T"='grey95', "N>G"='grey80', "N>C"='grey85', 
              "A>T"='orange2', "A>G"='orange3', "A>C"='orange4', "A>N"='orange1', 
              "T>A"='gold2', "T>G"='gold3', "T>C"='gold4', "T>N"='gold1',
              "G>A"='royalblue2', "G>T"='royalblue3', "G>C"='royalblue4', "G>N"='royalblue1', 
              "C>A"='steelblue2', "C>T"='steelblue3', "C>G"='steelblue4', "C>N"='steelblue1')
  
  posdata[, type := ordered(type, levels=levvec)]
  catlev <- rev(c("C:G in WRCY/RGYW", "C:G in other WRC/GYW", "C:G in SYC/GRS", "C:G in WRC/GYW and SYC/GRS", "other C:G", "A:T"))
  catcol <- rev(c("#440088", "#7700CC", "#AACCFF", "#AAAACC", "#0077FF", "#EE7700"))
  names(catcol) <- catlev
  # all(posdata$mutcat[posdata$mutated] %in% catlev)
  posdata[, mutcat := ordered(mutcat, levels=catlev)]

  
  # Calculate individual frequencies.
  maxdepth = max(posdata$depth)
  posdata[, freq := 0]
  posdata[(mutated), freq := count / depth]   # CANNOT use maxdepth as universal denominator because, in the merged tracks like Ramos, the library size of each fragment is different.
  
  # Remove mutations with high frequencies, as likely clonal SNPs.
  posdata <- posdata[!mutated | (mutated & freq <= allelecutoff), ]
  
  # Aggregate frequencies by position
  posdata[(mutated), aggrfreq := sum(.SD$freq), by=newpos]
  
  # Total mutations, so as to express the contingency tables as fraction of mutations.
  totmuts <- sum(posdata[(mutated), freq])
  
  # Mutation type frequencies
  posdata[mutated & !grepl('N', type), c('was', 'became') := list(substr(type, start=1, stop=1), 
                                                                  substr(type, start=3, stop=3))]
  # "From" as columns, "to" as rows
  contingency1 <- data.frame('to.from'=c('A', 'T', 'G', 'C', 'N'),
                             A=c(NA_real_,
                                 sum(posdata[mutated & was == 'A' & became == 'T', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'A' & became == 'G', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'A' & became == 'C', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'A' & became == 'N', freq]) / totmuts ),
                             T=c(sum(posdata[mutated & was == 'T' & became == 'A', freq]) / totmuts,
                                 NA_real_,
                                 sum(posdata[mutated & was == 'T' & became == 'G', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'T' & became == 'C', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'T' & became == 'N', freq]) / totmuts ),
                             G=c(sum(posdata[mutated & was == 'G' & became == 'A', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'G' & became == 'T', freq]) / totmuts,
                                 NA_real_,
                                 sum(posdata[mutated & was == 'G' & became == 'C', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'G' & became == 'N', freq]) / totmuts ),
                             C=c(sum(posdata[mutated & was == 'C' & became == 'A', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'C' & became == 'T', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'C' & became == 'G', freq]) / totmuts,
                                 NA_real_,
                                 sum(posdata[mutated & was == 'C' & became == 'N', freq]) / totmuts ),
                             N=c(sum(posdata[mutated & was == 'N' & became == 'A', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'N' & became == 'T', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'N' & became == 'G', freq]) / totmuts,
                                 sum(posdata[mutated & was == 'N' & became == 'C', freq]) / totmuts,
                                 NA_real_ ) )
  cont1 <- melt(as.data.table(contingency1), id.vars='to.from', variable.name='from', value.name='Frequency')
  setnames(cont1, c('to', 'from', 'Frequency'))
  cont1[, mutation := paste(from, to, sep='>')]
  cont1 <- cont1[!is.na(Frequency) & Frequency > 0,]
  cont1[, mutation := ordered(mutation, levels=levvec)]
  
  p1 <- ggplot(cont1, aes(x=mutation, y=Frequency*100, fill=mutation)) +
      geom_bar(stat='identity') +
      geom_hline(yintercept=0) +
      scale_y_continuous(trans='sqrt', expand=c(0,0), breaks=c(0,1,5,10,20,30,40,50,75,100)) +
      scale_fill_manual(values=colvec) +
      labs(x='', y='Mutations (%)') +
      theme_bw() +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
            legend.position='none',
            panel.border = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text=element_text(size=9),
            axis.title=element_text(size=9),
            legend.text=element_text(size=9),
            legend.title=element_blank() )
  p2 <- ggplot(cont1, aes(x=mutation, y=Frequency*100, fill=mutation)) +
      geom_bar(stat='identity', width=1, colour='white', size=0.2) +
      scale_y_sqrt() +
      scale_fill_manual(values=colvec) +
      coord_polar() +
      labs(x='', y='', title='Mutations (%)') +
      theme_bw() +
      theme(legend.position = 'none',
            axis.text=element_text(size=9),
            axis.title=element_text(size=9),
            legend.text=element_text(size=9),
            legend.title=element_blank() )
  p3 <- ggplot(cont1, aes(x='', y=Frequency*100, fill=mutation)) +
      geom_bar(stat='identity', position='stack', width=0.5, colour='white', size=0.2) +
      scale_fill_manual(values=colvec) +
      coord_polar('y', start=0) +
      theme_void() +
      theme(legend.position='right',
            legend.text=element_text(size=9),
            legend.title=element_blank() )
  
  # print( p1 + p2 + p3 )
  
  # Potential AIDER targets
  contingency2 <- data.frame(long_pattern_CG = sum(posdata[mutated & (sweetspotlong) & !sourspot, freq]) / totmuts,
                             short_pattern_CG = sum(posdata[mutated & (sweetspotshort) & !sourspot, freq]) / totmuts,
                             coldspot_CG = sum(posdata[mutated & sourspot & !(sweetspotlong | sweetspotshort), freq]) / totmuts,
                             ambiguous_CG = sum(posdata[mutated & sourspot & (sweetspotlong | sweetspotshort), freq]) / totmuts,
                             other_CG = sum(posdata[mutated & !(sweetspotlong | sweetspotshort) & !sourspot & grouptype=='G:C', freq]) / totmuts,
                             AT = sum(posdata[mutated & grouptype == 'A:T', freq]) / totmuts)
  cont2 <- as.data.table(t(contingency2))
  names(cont2) <- 'Frequency'
  cont2[, Category := names(contingency2)]
  cont2[Category=='long_pattern_CG', Category := 'C:G in WRCY/RGYW']
  cont2[Category=='short_pattern_CG', Category := 'C:G in other WRC/GYW']
  cont2[Category=='coldspot_CG', Category := 'C:G in SYC/GRS']
  cont2[Category=='other_CG', Category := 'other C:G']
  cont2[Category=='ambiguous_CG', Category := 'C:G in WRC/GYW and SYC/GRS']
  cont2[Category=='AT', Category := 'A:T'] 
  cont2[, Category := ordered(Category, levels=catlev)]
  
  # Category clump
  clump <-  unique(posdata[(mutated), .(pos, aggrfreq, mutcat)])
  setorder(clump, -aggrfreq, na.last = TRUE)
  clump[, pos := ordered(pos, levels=pos)]
  
  p4 <- ggplot(cont2, aes(x='', y=Frequency*100, fill=Category) ) +
    geom_bar(stat='identity', position='stack', colour='white') +
    # geom_hline(yintercept = 0) +
    scale_fill_manual(values=catcol) +
    labs(x='', y='', title='Mutations (%)') +
    theme_bw() +
    theme(legend.position='right',
          axis.ticks.length.x = unit(0, 'cm'),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.text=element_text(size=9),
          axis.title=element_text(size=9),
          legend.text=element_text(size=9),
          legend.title=element_blank())
  p5 <- ggplot(cont2, aes(x='', y=Frequency*100, fill=Category) ) +
    geom_bar(stat='identity', position='stack', colour='white') +
    # geom_hline(yintercept = 0) +
    scale_fill_manual(values=catcol) +
    coord_polar('y', start=0) +
    labs(x='', y='', title='Mutations (%)') +
    theme_void() +
    theme(legend.position='none',
          legend.text=element_text(size=9),
          legend.title=element_blank() )
  
  p6 <- ggplot(clump, aes(x=pos, y=aggrfreq, fill=mutcat)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=catcol) +
    labs(x='Positions', y='Mutation Frequency') +
    theme_minimal() +
    theme(legend.position='right',
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title=element_text(size=9),
          axis.ticks.length.x = unit(0, 'cm'))
  
  # print( p4 + p5 + p6)
  
  print(p1 + p2 + p4)
  print(p6)
  
  # Write out mutation type counts
  cat(c(basename(sf), "\n\n"), file=mtxf, append=TRUE)
  write.table(contingency1, file=mtxf, sep="\t", row.names=FALSE, quote=FALSE, append=TRUE)
  cat("\n", file=mtxf, append=TRUE)
  write.table(contingency2, file=mtxf, sep="\t", row.names=FALSE, quote=FALSE, append=TRUE)
  cat("\n", file=mtxf, append=TRUE)
  
  # Top 5% highest peaks
  posdata[(mutated), flag := aggrfreq >= quantile(aggrfreq, 0.95)]
  
  # Fill in gaps in the coordinates.
  setorder(posdata, newpos)
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
                                  newpos = (posdata[missing[i, from], newpos] + 1):(posdata[missing[i, from], newpos] + n),
                                  seq = rep(posdata[missing[i, from], seq], n),
                                  mutated = FALSE,
                                  grouptype = NA_character_,
                                  sweetspotlong = NA,
                                  sweetspotshort = NA,
                                  sourspot = NA,
                                  mutcat = NA_character_,
                                  freq = 0,
                                  aggrfreq = NA_real_,
                                  was = NA_character_,
                                  became = NA_character_,
                                  flag = NA) )
    }
  }
  

  # Plot mutation frequencies by position
  
  DATA1 <- unique( posdata[(mutated & aggrfreq >= 0), .(seq, newpos, aggrfreq, mutcat, depth)] )
  # Shift data coordinates so the x axis starts at 1
  DATA1[, newpos := newpos - showfrom + 1]
  
  minfreq <- min(DATA1$aggrfreq, na.rm=TRUE)
  maxfreq <- ifelse(any(DATA1$aggrfreq > 0),
                    max(DATA1[(aggrfreq > 0), aggrfreq], na.rm=TRUE),  # Use the largest positive value, if any
                    abs(max(DATA1[(aggrfreq >= 0), aggrfreq], na.rm=TRUE)) ) # use symmetric to the smallest negative value
  
  # Base
  pp <- ggplot(DATA1,
               aes(x = newpos, y = aggrfreq, fill = mutcat)) +
    facet_grid(seq ~ .) +
    scale_x_continuous(expand=c(0,1), breaks=seq(20, showto - showfrom + 1, 20)) +
    labs(x="Position") +
    theme(legend.position = 'bottom',
          legend.title = element_blank(),
          axis.text=element_text(size=9),
          axis.title=element_text(size=9),
          panel.background = element_rect(fill="white"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  # Annotation
  covscale2 <- 1.1
  pplain <- pp
  if (nrow(cdrs) > 0) {
    for (i in 1:nrow(cdrs)) {
      pplain <- pplain +
        geom_rect(xmin=cdrs[i,]$xleft, xmax=cdrs[i,]$xright, ymin=-0.2 * maxfreq, ymax=covscale2 * maxfreq, fill='grey95')
    }
  }
  if (nrow(warmspotshort) > 0) {
    pplain <- pplain +
      geom_segment(data=warmspotshort,
                   aes(x=xleft, xend=xright, y= -0.1 * maxfreq, yend= -0.1 * maxfreq), inherit.aes=FALSE,
                   colour='grey70', size=3)
  }
  if (nrow(hotspotshort) > 0) {
    pplain <- pplain +
      geom_segment(data=hotspotshort,
                   aes(x=xleft, xend=xright, y= -0.1 * maxfreq, yend= -0.1 * maxfreq), inherit.aes=FALSE,
                   colour='red1', size=3)
  }
  if (nrow(warmspotslong) > 0) {
    pplain <- pplain +
      geom_segment(data=warmspotslong,
                   aes(x=xleft, xend=xright, y= -0.05 * maxfreq, yend= -0.05 * maxfreq), inherit.aes=FALSE,
                   colour='grey50', size=3)
  }
  if (nrow(hotspotslong) > 0) {
    pplain <- pplain +
      geom_segment(data=hotspotslong,
                   aes(x=xleft, xend=xright, y= -0.05 * maxfreq, yend= -0.05 * maxfreq), inherit.aes=FALSE,
                   colour='red3', size=3)
  }
  if (nrow(coldspots) > 0) {
    pplain <- pplain +
      geom_segment(data=coldspots,
                   aes(x=xleft, xend=xright, y= -0.15 * maxfreq, yend= -0.15 * maxfreq), inherit.aes=FALSE,
                   colour='black', size=3)
  }
  # Content
  pplain3 <- pplain +
    geom_bar(data=DATA1, stat = 'identity', width=1, colour='white', size=rel(0.2)) +
    geom_hline(yintercept=0, size=0.6) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=catcol, na.value="transparent") +
    coord_cartesian(xlim=c(1, showto - showfrom + 1),
                    ylim=c(-0.2 * maxfreq, maxfreq * covscale2)) +
    labs(y='Frequency', title=basename(sf)) +
    theme(legend.position = 'none')
  
  # Also plot out the coverage profile. This is the subscript for the frequencies. Large coverage can cause misleading frequencies.
  # Substracted tracks, keeping the lowest coverage available, maybe especially prone to dips.
  md <- max(DATA1$depth)
  pcov <- pplain +
    geom_line(aes(x=newpos, y=depth), inherit.aes = FALSE) +
    scale_y_continuous(limits=c(0, 10 ^ floor(log10(md)) + 10 ^ floor(log10(md)) * ceiling(md / 10^ floor(log10(md))) ),
                       expand=c(0,0)) +
    labs(y='Read coverage above Base Quality cutoff', title=basename(sf))
  
  print(pcov)
  print( pplain3 )
  
  # Create bedGraph of mutations
  # Use the normal mapped coordinates for the bedGraph, not the shifted ones used for peak labelling in the plots.
  bG <- unique(posdata[(mutated), .(chr, pos-1, pos, aggrfreq)])
  if (dim(bG)[1] != 0) {
    bf <- file.path(outdir, basename(sub('stats', paste0('cov', maxdepth, '.bedGraph'), sf)))
    cat(paste0("track type=bedGraph name=", basename(sf), "\n"), file=bf)
    fwrite(bG, file=bf, col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE, append=TRUE)
  } else {
    message("No mutations to put in the bedGraph.")
  }
  
}

dev.off()

