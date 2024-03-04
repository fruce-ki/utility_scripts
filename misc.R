
# RMarkdown header template
###########################

---
title: "Report"
author: "Kimon Froussios"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
    code_folding: hide
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: false
editor_options:
  chunk_output_type: console
params:
  foo: "bar"
---



# Combinations
G <- as.data.table(t( combn(samples, 2) ))        # unique pairs
G <- as.data.table(expand.grid(samples, samples)) # includes mirrored pairs and self-pairs

### Colour Palettes
###################
{
  library(ggplot2)

  # convert list of RGB triplets. to hex
  rgb2hex <- function(RGBtuples){ rgb( t(as.data.frame(RGBtuples)), maxColorValue = 255) }

  # convert colour name vector to hex
  colour2hex <- function(namevector) { rgb( t(as.data.frame( lapply(namevector, function(x){ col2rgb(x)/255 } ) )) ) }

  # preview colour vector
  showpalette <- function(p) {
    p <- factor(p, ordered=TRUE, levels=p)
    ggplot(data.frame(x=p, y=1), aes(p, y, fill=p)) +
      geom_bar(stat='identity', colour='transparent') +
      scale_fill_identity() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
            legend.position = 'none')
  }

  # colourblind palette
  palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  showpalette(palette_OkabeIto)

  install.packages("rcartocolor")
  library(rcartocolor)
  scales::show_col(carto_pal(12, "Safe"))

  ramp_colours <- c(
    "#bb0000", "#eedd00", "#008800", "#00dddd", "#0000ff", "#ff00ff"
  )
  showpalette(ramp_colours)

  mycolours <- c(
    "#880000", "#bb2200", "#ff0000",
    "#bb6600", "#ff8833", "#eeaa00", "#ffcc00",
    "#0000ff", "#0066ff", "#0099ff", "#00ccff",  "#00ffff",
    "#5500aa", "#8800ff", "#8866ff", "#ff00dd", "#ff77ff",  "#eeaaff", "#ffccff",

    "#224400", "#226600", "#229900", "#22cc00", "#22ff00", "#88cc00",
    "#8899aa",
    "#0033aa", "#0066aa", "#0099aa", "#00ccaa",
    "#886600", "#bb9900", "#889900",
    "#8800aa", "#8866aa", "#8899ff", "#88ccff",
    "#ff99bb", "#eeaabb", "#eeeeaa"
  )

  categorical_colours <- c(
    "#bb2200", "#0066ff", "#8800ff", "#ff8833", "#229900",
    "#880000", "#0000ff", "#5500aa", "#bb6600", "#224400",
    "#eeaabb", "#00ccff", "#ff77ff", "#ffcc00", "#22cc00",
    "#8899aa",
    "#ff00dd", "#0099ff", "#8866ff", "#22ff00", "#886600",
    "#889900", "#88ccff", "#eeeeaa", "#0066aa", "#0099aa"
  )

  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )

  # install.packages("Polychrome")
  library(Polychrome)

  # build-in color palette
  Glasbey = glasbey.colors(32)
  swatch(Glasbey)
}


### Rounding to any precision
#############################
mround <- function(x, base){
  base*round(x/base)
}


### Rolling slice
#################
# Returns a list of index vectors for x.
rollslice <- function(x, n) {
  # x <- c(1,2,3,4,5,6,7,8,9,0)
  # n <- 4
  l <- length(x)
  s <- as.integer(l/n)  # number of slices
  if (s * n < l) {
    s <- s + 1
  }

  lapply(1:s, function(y){
    # y <- 2    # 1 -> 1:4,  2 -> 5:8,  3 -> 9:10
    ((y-1)*n + 1):min(l, y*n)
  })
}


### Geometric Mean
###################
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


### Correlations
################
my_pairwise_internal_corels <- function(mat, samples, method = "pearson", rds=NULL, txs=3, minMean=0, minSingle=0, groups=NULL, loopVal="default") {
  # mat <- log10(counts[topgenes,]);

  # Filter
  if (minMean != 0 | minSingle != 0) {
    mat <- mat[rowSums(mat >= minSingle) >= 1 | rowMeans(mat) >= minMean, ]
  } else {
    mat <- mat[rowSums(mat) > 0, ]
  }

  # Correlations
  cormat <- cor(mat, method=method)

  # Cluster
  hcfit <- hclust(dist(scale(cormat, center=TRUE)))
  rn <- rownames(cormat)

  # Make dendrogram. https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
  dend <- as.dendrogram(hcfit)
  dend_data <- dendro_data(dend)
  # Setup the data, so that the axes are exchanged, instead of using coord_flip()
  segment_data <- with(
    segment(dend_data),
    data.frame(x = y, y = x, xend = yend, yend = xend))
  # Use the dendrogram label data to position the sample labels
  sample_pos_table <- with(
    dend_data$labels,
    data.frame(y_center = x, sample = as.character(label), height = 1))
  # Limits for the vertical axes
  sample_axis_limits <- with(
    sample_pos_table,
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
  ) + 0.1 * c(-1, 1) # extra spacing: 0.1
  # Dendrogram plot
  pd <- ggplot(segment_data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_reverse(expand = c(0, 0.5),
                    position = "top") +
    scale_y_continuous(position = "right",
                       breaks = sample_pos_table$y_center,
                       labels = sample_pos_table$sample,
                       limits = sample_axis_limits,
                       expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())


  # Create duplicates for different plot styles.
  cormat <- cormat[samples, samples]                     # In the supplied order.
  cormat_t <- cormat                                     # Duplicate in which to delete below the diagonal.
  cormat_c <- cormat[rn[hcfit$order], rn[hcfit$order]]   # Duplicate in clustered order.
  cormat_ct <- cormat_c                                  # Duplicate in which to delete below the diagonal.
  # Delete below the diagonal, for the triangles.
  for (r in 1:nrow(cormat_t)) {
    for (c in 1:ncol(cormat_t)) {
      if (c <= r) {                # For non-clustered, also delete the diagonal.
        cormat_t[r, c] <- NA_real_
      }
    }
  }
  for (r in 1:nrow(cormat_ct)) {
    for (c in 1:ncol(cormat_ct)) {
      if (c < r) {                 # For clustered keep the diagonal, so the dendrogram lines up.
        cormat_ct[r, c] <- NA_real_
      }
    }
  }
  cormat_ctv <- cormat_ct                               # Duplicate in which to also delete the diagonal, for the value labels.
  for (r in 1:nrow(cormat_ctv)) {
    for (c in 1:ncol(cormat_ctv)) {
      if (c == r) {
        cormat_ctv[r, c] <- NA_real_
      }
    }
  }

  # Restructure for plotting.
  restruct <- function(cm){
    rn <- rownames(cm)
    cm <- as.data.table(cm)
    cm[, observation1 := factor(rn, ordered=TRUE, levels=rn)]
    cm <- melt(cm, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
    cm[, observation2 := factor(observation2, ordered=TRUE, levels=rn)]
    # cm <- merge(cormat, sample_pos_table, by.x="observation2", by.y="sample", all.x=TRUE)
    cm
  }

  cormat <- restruct(cormat)
  cormat_t <- restruct(cormat_t)
  cormat_c <- restruct(cormat_c)
  cormat_ct <- restruct(cormat_ct)
  cormat_ctv <- restruct(cormat_ctv)
  cormat_t <- cormat_t[!is.na(Correlation)]
  cormat_ct <- cormat_ct[!is.na(Correlation)]
  cormat_ctv <- cormat_ctv[!is.na(Correlation)]


  # Text colour switch for the dynamic range
  m <- min(cormat4$Correlation, na.rm=TRUE)
  M <- max(cormat4$Correlation, na.rm=TRUE)
  colourswitch <- c( m + 0.49 * (M-m),  m + 0.51 * (M-m) )


  # Square. Custom order. No values. Full range.
  p_sonf <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "grey50" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )

  # Square. Custom order. No values. Dynamic range.
  p_sond <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "grey50" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )


  # Square. Custom order. With values. Full range.
  p_sovf <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat_t, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "grey50" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="grey50", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )

  # Square. Custom order. With values. Dynamic range.
  p_sovd <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat_t, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "grey50" ) +
    scale_colour_manual(values=c("black", "white"), na.value="grey50", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )


  # Square. Clustered order. No values. Full range.
  p_scnf <- ggplot(cormat_c, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "grey50" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )

  # Square. Clustered order. No values. Dynamic range.
  p_scnd <- ggplot(cormat_c, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "grey50" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )


  # Square. Clustered order. With values. Full range.
  p_scvf <- ggplot(cormat_c, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat_ctv, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "grey50" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="grey50", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )

  # Square. Clustered order. With values. Dynamic range.
  p_scvd <- ggplot(cormat_c, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat_ctv, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "grey50" ) +
    scale_colour_manual(values=c("black", "white"), na.value="grey50", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )


  # Triangle. Custom order. With values. Full range.
  p_tovf <- ggplot(cormat_t, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "grey50" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="grey50", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )

  # Triangle. Custom order. With values. Dynamic range.
  p_tovd <- ggplot(cormat_t, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "grey50" ) +
    scale_colour_manual(values=c("black", "white"), na.value="grey50", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )


  # Triangle. Clustered order. with values. Full range.
  p_tcvf <- ggplot(cormat_ct, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat_ctv, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "grey50" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="grey50", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )

  # Triangle. Clustered order. with values. Dynamic range.
  p_tcvd <- ggplot(cormat_ct, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat_ctv, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value ="grey50" ) +
    scale_colour_manual(values=c("black", "white"), na.value="grey50", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )


  # If exact pairs don't matter.
  cormat_x <- merge(cormat_t, groups, by = c('observation1', 'observation2'), all.x = TRUE, all.y = FALSE)

  p_bee <- ggplot(cormat_x, aes(x = Correlation, y = group, colour = group)) +
    # geom_violin()+ #(draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(outlier.alpha = 0, width = 0.5, fill = 'grey95', colour = 'black') +
    # geom_jitter(width = 0, height = 0.2) +
    geom_beeswarm(groupOnX = FALSE) +
    # scale_x_continuous(limits = c(-1, 1)) +
    scale_colour_brewer(palette = 'Set1') +
    labs(title = paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Summary"), x = "Pairwise Correlation", y = "Pairs") +
    theme(legend.position = 'none')

  out <- list(corr=dcast(cormat2, observation1 ~ observation2, value.var = "Correlation"),
              sonf=p_sonf, #sond=p_sond,
              sovf=p_sovf, #sovd=p_sovd,
              scnf=pd + p_scnf + plot_layout(ncol=2, widths=c(1,4)), #scnd=pd + p_scnd + plot_layout(ncol=2, widths=c(1,4)),
              scvf=pd + p_scvf + plot_layout(ncol=2, widths=c(1,4)), #scvd=pd + p_scvd + plot_layout(ncol=2, widths=c(1,4)),
              tovf=p_tovf, #tovd=p_tovd,
              tcvf=pd + p_tcvf + plot_layout(ncol=2, widths=c(1,4)), #tcvd=pd + p_tcvd + plot_layout(ncol=2, widths=c(1,4)),
              bee = p_bee
  )

  if (!is.null(rds) && length(rds) > 0) {
    saveRDS(out, file = rds)
  }

  return(out)
}




### Bisque
##########
parse_bisque <- function(input = params$path) {
  stopifnot(require(data.table))

  cat(input, "\n")
  DT <- data.table::fread(input)
  setnames(DT, c('Sample', names(DT)[2:length(DT)]))
  setkey(DT, Sample)

  input <- sub('.tsv$', '', basename(input))
  input <- sub('_compositions', '', input)
  input <- sub('Bisque_', '', input)
  input <- strsplit(input, '_')[[1]]
  DT[, study := input[1]]
  DT <- DT[!grepl('^ *$', Sample), ]
  DT <- melt(DT, id.vars = c('Sample', 'study'), variable.name = 'category', value.name = 'fraction')
  if (length(input) == 2)
    DT[, category := paste(input[2], category, sep = '__')]

  return(DT)
}


### TPM/RPM
###########
normscale <- function(counts, featLens = NULL, specnorm = NULL){
  # counts <- M
  TPMs <- NULL
  if (!is.null(featLens)) {
    stopifnot(all(is.finite(featLens)))
    stopifnot(length(featLens) == nrow(counts))
    # Scale by feature size
    message("Scale to TPM.")
    TPMs <- sweep(counts, 1, featLens, `/`)
  } else {
    message("Scale to RPM.")
    TPMs <- counts
  }
  # Scale by sequencing depth
  if (!is.null(specnorm)) {
    message("Scale with special exclusions.")
    colsums <- colSums(TPMs[!grepl(specnorm, rownames(TPMs), ignore.case=TRUE, perl=TRUE), ], na.rm=TRUE)
  } else {
    colsums <- colSums(TPMs, na.rm=TRUE)
  }
  TPMs <- sweep(TPMs, 2, colsums, `/`) * 1e6

  return(TPMs)
}


### Gene length
###############
genlen_from_coord <- function(coord, right_closed = 1) {
  # coord <- c('chr1:11869-14409:+', 'chrF:11111-22222:-')
  a <- vapply(strsplit(coord, ':'), function(x) { x[2] }, character(1))
  b <- vapply(strsplit(a, '-'), function(x) { as.integer(x[2]) - as.integer(x[1]) + right_closed }, numeric(1))

  return(b)
}


### GGpairs for continuous variables within a dataframe.
########################################################
{
  # Coefficients
  my_cont_ggpairs_upperFn <- function(data, mapping, info="corr", ...) {
    x <- sub('~', '', as.character(mapping["x"]))
    y <- sub('~', '', as.character(mapping["y"]))
    bingo <- NA
    if (info == "corr"){
      cm <- cor(data)
      bingo <- round(cm[y, x], 2)
    } else if (info == "lmrsq") {
      l <- lm(data[[y]] ~ data[[x]])
      bingo <- round(summary(l)$adj.r.square, 2)
    } else {
      stop()
    }
    cpal <- colorRampPalette(c("darkred", "darkorange", "gold", "steelblue", "darkblue"))(201)
    p <- ggplot() +
      annotate("text", x=0.5, y=0.5, label=bingo, colour=cpal[bingo * 100 + 101], ...) +
      coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
      theme(panel.background = element_rect(fill="white"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())
    return(p)
  }

  # Scatters
  # *trans are normal ggplot2 transformation strings like "log10"
  my_cont_ggpairs_lowerFn <- function(data, mapping, method = "lm", xtrans=NULL, ytrans=NULL, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_abline(intercept = 0, slope=1) +
      geom_point(colour = "blue", shape=16, size=rel(0.8), alpha=0.2) +
      geom_smooth(method = method, color = "red", ...) +
      theme(panel.background = element_rect(fill="white"))
    if (!is.null(xtrans)) {
      p <- p + scale_x_continuous(trans=xtrans)
    }
    if (!is.null(ytrans)) {
      p <- p + scale_y_continuous(trans=ytrans)
    }
    return(p)
  }

  # Distributions
  my_cont_ggpairs_diagFn <- function(data, mapping, xtrans=NULL, ...) {
    x <- sub('~', '', as.character(mapping["x"]))
    # v <- data[, x, with=FALSE]
    # names(v) <- 'x'
    p <- ggplot(data, mapping) +
      geom_density(...) +
      theme(panel.background = element_rect(fill="grey95"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())
    if (!is.null(xtrans)) {
      p <- p + scale_x_continuous(trans=xtrans)
    }
    return(p)
  }

  # Put it all together (because I'll forget the syntax)
  my_cont_ggpairs <- function(df, method="lm", info="lmrsq", xtrans=NULL, ytrans=NULL) {
    p <- ggpairs(df,
          lower = list(continuous = wrap(my_cont_ggpairs_lowerFn, method=method, xtrans=xtrans, ytrans=ytrans)),
          # lower = list(continuous = "na"),
          diag = list(continuous = wrap(my_cont_ggpairs_diagFn, xtrans=xtrans)),
          # diag = list(continuous = "naDiag"),
          upper = list(continuous = wrap(my_cont_ggpairs_upperFn, info=info, size=5))
      ) +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
    return(p)
  }

}


