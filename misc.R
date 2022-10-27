


# RMarkdown header template
###########################

# title: "Report"
# author: "Kimon Froussios"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document:
#     code_folding: hide
#     toc: true
#     toc_float: true
#     toc_depth: 4
#     collapsed: false
# params:


### Colour Palettes
###################

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
showpalette( c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") )

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
showpalette(mycolours)

mycolours <- c(
  "#bb2200", "#0066ff", "#8800ff", "#ff8833", "#229900", 
  "#880000", "#0000ff", "#5500aa", "#bb6600", "#224400",
  "#eeaabb", "#00ccff", "#ff77ff", "#ffcc00", "#22cc00",
  "#8899aa",
  "#ff00dd", "#0099ff", "#8866ff", "#22ff00", "#886600", 
  "#889900", "#88ccff", "#eeeeaa", "#0066aa", "#0099aa"
)
showpalette(mycolours)
length(mycolours)


### Rounding to nearest non-10
##############################
mround <- function(x,base){
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



cormat <- matrix(rep(seq(-0.9, 1, 0.1),5), ncol = 10)
samples <- c('a','b','c','d','e','f','g','h','i','j')
colnames(cormat) <- samples
rownames(cormat) <- samples

### Correlations within the dataframe.
######################################
# colnames and rownames yes, non-numeric columns no.
# Requires a vector of sample names for the non-clustered plots.
# Correlations
my_pairwise_internal_corels <- function(mat, samples, method = "pearson", rds=NULL, txs=3, minMean=0, minSingle=0, sizefactor=1) {
  # mat <- counts; samples <- covars$Samp
  
  # if (method == "pearson") 
  libsizes <- colSums(mat)
  # Filter
  if (minMean != 0 | minSingle != 0) {
    mat <- mat[rowSums(mat >= minSingle) >= 1 | rowMeans(mat) >= minMean, ]
  } else {
    mat <- mat[rowSums(mat) > 0, ]
  }
  # Scale
  # if (method == "pearson") 
  mat <- sweep(mat, 2, libsizes, `/`) * sizefactor
  
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
    data.frame(y_center = x, Sample = as.character(label), height = 1))
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
                       labels = sample_pos_table$Sample, 
                       limits = sample_axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = NULL, y = NULL) +
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  
  # Create duplicates for different plot styles
  cormat <- cormat[samples, samples]                    # Supplied order
  cormat2 <- cormat                                     # Duplicate in which to delete below the diagonal.
  cormat3 <- cormat[rn[hcfit$order], rn[hcfit$order]]   # Duplicate in clustered order.
  cormat4 <- cormat3                                    # Duplicate in clustered order in which to delete below the diagonal.
  # Delete below diagonal half for the numeric labels.
  for (r in 1:nrow(cormat2)) {
    for (c in 1:ncol(cormat2)) {
      if (c <= r) {
        cormat2[r, c] <- NA_real_
      }
    }
  }
  for (r in 1:nrow(cormat4)) {
    for (c in 1:ncol(cormat4)) {
      if (c <= r) {
        cormat4[r, c] <- NA_real_
      }
    }
  }
  
  # Restructure for plotting.
  rn <- rownames(cormat)
  cormat <- as.data.table(cormat)
  cormat[, observation1 := factor(rn, ordered=TRUE, levels=rn)]
  cormat <- melt(cormat, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat[, observation2 := factor(observation2, ordered=TRUE, levels=rn)]
  # cormat <- merge(cormat, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
  
  rn2 <- rownames(cormat2)
  cormat2 <- as.data.table(cormat2)
  cormat2[, observation1 := factor(rn2, ordered=TRUE, levels=rn2)]
  cormat2 <- melt(cormat2, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat2[, observation2 := factor(observation2, ordered=TRUE, levels=rn2)]
  cormat2 <- cormat2[!is.na(Correlation)]
  # cormat2 <- merge(cormat2, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
  
  rn3 <- rownames(cormat3)
  cormat3 <- as.data.table(cormat3)
  cormat3[, observation1 := factor(rn3, ordered=TRUE, levels=rn3)]
  cormat3 <- melt(cormat3, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat3[, observation2 := factor(observation2, ordered=TRUE, levels=rn3)]
  # cormat3 <- merge(cormat3, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
  
  rn4 <- rownames(cormat4)
  cormat4 <- as.data.table(cormat4)
  cormat4[, observation1 := factor(rn3, ordered=TRUE, levels=rn4)]
  cormat4 <- melt(cormat4, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat4[, observation2 := factor(observation2, ordered=TRUE, levels=rn4)]
  cormat4 <- cormat4[!is.na(Correlation)]
  # cormat4 <- merge(cormat4, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
  
  # Text colour switch for the dynamic range
  m <- min(cormat4$Correlation, na.rm=TRUE)
  M <- max(cormat4$Correlation, na.rm=TRUE)
  colourswitch <- c( m + 0.49 * (M-m),  m + 0.51 * (M-m) )
  
  
  # # Square. Custom order. No values. Full range.
  # p1 <- ggplot(cormat, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
  #   scale_x_discrete(position = "top") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
  #   theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  # 
  # # Square. Custom order. No values. Dynamic range.
  # p1a <- ggplot(cormat, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
  #   scale_x_discrete(position = "top") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
  #   theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  # 
  # # Triangle. Custom order. With values. Full range.
  # p2 <- ggplot(cormat2, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
  #   scale_x_discrete(position = "top") +
  #   scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
  #   scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
  #   theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  # 
  # # Triangle. Custom order. With values. Dynamic range.
  # p2a <- ggplot(cormat2, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
  #   scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
  #   scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
  #   scale_x_discrete(position = "top") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
  #   theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  
  # Square. Custom order. With values. Full range.
  p12 <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat2, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Custom order. With values. Dyhamic range.
  p12a <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat2, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
    scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  
  # # Square. Clustered order. No values. Full range.
  # p3 <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
  #   scale_x_discrete(position = "top") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
  #   theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  # 
  # # Square. Clustered order. No values. Dyhamic range.
  # p3a <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
  #   scale_x_discrete(position = "top") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
  #   theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  # 
  # # Triangle. Clustered order. with values. Full range.
  # p4 <- ggplot(cormat4, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
  #   scale_x_discrete(position = "top") +
  #   scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
  #   scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
  #   theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  # 
  # # Triangle. Clustered order. with values. Dyhamic range.
  # p4a <- ggplot(cormat4, aes(x=observation1, y=observation2)) +
  #   geom_tile(aes(fill=Correlation)) +
  #   geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
  #   scale_x_discrete(position = "top") +
  #   scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
  #   scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
  #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
  #   theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
  #         panel.grid = element_blank() )
  
  # Square. Clustered order. With values. Full range.
  p34 <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat4, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Clustered order. With values. Dyhamic range.
  p34a <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat4, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(aspect.ratio=1, axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  out <- list(corr=dcast(cormat2, observation1 ~ observation2, value.var = "Correlation"),
              # sfrnc=p1, tfrnc=p2, sdrnc=p1a, tdrnc=p2a,
              frnc=p12, drnc=p12a,
              # sfrc=pd + p3, tfrc=p4, sdrc=pd + p3a, tdrc=p4a,
              frc=pd + p34 + plot_layout(ncol=2, widths=c(1,4)), 
              drc=pd + p34a + plot_layout(ncol=2, widths=c(1,4)))
  
  if (!is.null(rds)) {
    saveRDS(out, file = rds)
  }
  
  return(out)
}


### Correlations in a sparse matrix.
####################################
sparse.cor <- function(x){
  # https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
  
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element
  
  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans
  
  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
                crossprod(t(Ex))*(n-length(ii))
  )/(n-1)
  sdvec <- sqrt(diag(covmat))
  covmat/crossprod(t(sdvec))
}


# PCA
# (scaling and centering of genes to one another, NOT samples to one another)
pca_plotter <- function(pc, loads, highloads, ig, pcaimp, pcx, pcy){  ## helper function
  # pcx <- 1; pcy <- 2
  mycolours <- c(
    "#dd2200", "#0066ff", "#8800ff", "#ff8833", "#229900", 
    "#880000", "#0000ff", "#5500aa", "#bb6600", "#224400",
    "#eeaabb", "#00ccff", "#ff77ff", "#ffcc00", "#22cc00",
    "#8899aa",
    "#ff00dd", "#0099ff", "#8866ff", "#22ff00", "#886600", 
    "#889900", "#88ccff", "#eeeeaa", "#0066aa", "#0099aa"
  )
  pcxs <- paste0("PC", pcx)
  pcys <- paste0("PC", pcy)
  subload <- loads[rowID %in% unique(c(head(highloads[[pcx]], 5), head(highloads[[pcy]], 5))), c("rowID", pcxs, pcys), with=FALSE]
  lmt <- max( c(abs(pc[[pcxs]]), abs(pc[[pcys]])) )
  
  pc12 <- wrap_plots(lapply(ig, function(varname) {
    return(
      ggplot(pc, aes_string(x=pcxs, y=pcys, label="Sample", fill=varname, colour=varname)) +
        geom_vline(xintercept=0, linetype='dashed', colour='grey70') +
        geom_hline(yintercept=0, linetype='dashed', colour='grey70') +
        geom_mark_hull(alpha=0.1, colour='transparent') +
        geom_point() +
        # geom_text_repel(data=subload, inherit.aes=FALSE, aes_string(label="rowID", x=pcxs, y=pcys), size=2) +
        coord_fixed(ratio=1, xlim=c(-lmt, lmt), ylim=c(-lmt, lmt) ) +
        scale_colour_manual(values=mycolours) +
        scale_fill_manual(values=mycolours) +
        labs(x=paste0(pcxs, " (", round(pcaimp[PC==pcx, Explained], 1), "%)"), 
             y=paste0(pcys, " (", round(pcaimp[PC==pcy, Explained], 1), "%)")) +
        theme(panel.grid=element_blank())
    )}), ncol = 2)
  
  pc12v <- wrap_plots(lapply(ig, function(varname) {
    return(
      ggplot(pc, aes_string(x=pcxs, y=pcys, label="Sample", fill=varname, colour=varname)) +
        geom_vline(xintercept=0, linetype='dashed', colour='grey70') +
        geom_hline(yintercept=0, linetype='dashed', colour='grey70') +
        # geom_mark_hull(alpha=0.1, colour='transparent') +
        geom_segment(data=subload, inherit.aes=FALSE, colour='black', arrow=arrow(length=unit(0.1, 'cm')), alpha=0.2,
                     x=0, y=0, aes_string(xend=pcxs, yend=pcys)) +
        geom_point() +
        geom_text_repel(data=subload, inherit.aes=FALSE, aes_string(label="rowID", x=pcxs, y=pcys), size=2) +
        coord_fixed(ratio=1, xlim=c(-lmt, lmt), ylim=c(-lmt, lmt) ) +
        scale_colour_manual(values=mycolours) +
        scale_fill_manual(values=mycolours) +
        labs(x=paste0(pcxs, " (", round(pcaimp[PC==pcx, Explained], 1), "%)"), 
             y=paste0(pcys, " (", round(pcaimp[PC==pcy, Explained], 1), "%)")) +
        theme(panel.grid=element_blank())
    )}), ncol = 2)
  
  return(list(pc12, pc12v))
}

do_pca <- function(countsmat, covars, scale = TRUE, center = TRUE, rds=NULL, topgenes=NULL, minMean=0, minSingle=0, ntop=10) {
  # countsmat <- counts; covars <- covars
  
  nvars <- nrow(countsmat)
  countsmat <- countsmat[rowSums(countsmat) > 0, ]
  
  # Mean and Standard Deviation, standardized counts, mean and stdev of standardized counts.
  genevar <- data.table(name = rownames(countsmat),
                        Mean = rowMeans(countsmat),
                        StDev = rowSds(countsmat),
                        singles = rowSums(countsmat > minSingle))
  normcounts <- (countsmat - genevar$Mean) / genevar$StDev
  genevar [, zMean := rowMeans(normcounts)]
  # genevar [, zStDev := rowSds(normcounts)]    # always 1 by definition
  
  # Select features
  if (is.null(topgenes)) {
    message(paste("Using variable features that exceed either", minMean, "mean count or", minSingle, "count in any single sample."))
    topgenes = genevar[StDev > 0 & (Mean >= minMean | singles >= 1), name]  # all the variable genes, above min count level
  } else {
    message("Using variable features among the manually provided list.")
    topgenes = genevar[name %in% topgenes & StDev > 0, name] # manually provided features, as long as they are variable
  }
  subcmat <- countsmat[topgenes, ]
  
  # Rotate counts so genes are variables and samples are observations
  subcmat <- t(subcmat)
  message(paste("Number of variable features available:", nrow(genevar[StDev>0,])))
  message(paste("Number of variable features used:", length(topgenes)))
  
  # PCA
  pca <- prcomp(subcmat, center = center, scale = scale)
  
  srn <- sqrt(nrow(pca$x) - 1)
  pc <- sweep(pca$x, 2, 1 / (pca$sdev * srn), FUN = '*')   
  dirs <- as.data.table(pca$rotation)
  
  pc <- cbind(pc, data.frame(Sample = rownames(pc)))
  pc <- as.data.frame(merge(pc, covars, by="Sample", all = TRUE))
  npc <- sum(pca$sdev > 1)
  
  # Screeplot.
  
  pcaimp <- data.table(PC = 1:length(colnames(pca$x)),
                       Explained = summary(pca)$importance['Proportion of Variance', ] * 100,
                       Cumulative = summary(pca)$importance['Cumulative Proportion', ] * 100)
  pcaimp <- pcaimp[1:npc, ]
  topnpc <- nrow(pcaimp[Explained >= 1])
  
  pimp <- ggplot(pcaimp) +
    geom_line(aes(x=PC, y=Cumulative), colour='dodgerblue') +
    geom_bar(aes(x=PC, y=Explained, fill=Explained >= 1),
             stat='identity', colour='transparent', alpha=0.3) +
    geom_text(aes(x=PC, y=0, label=paste0(round(Explained, 1),"%")), 
              angle=90, hjust=0, vjust=0, size=rel(3)) +
    geom_text(data=pcaimp[2:npc,],
              aes(x=PC, y=0.95*Cumulative, label=paste0(round(Cumulative, 1),"%")), 
              angle=90, vjust=0, hjust=1, size=rel(3)) +
    scale_fill_manual(values=c("grey25", "grey75"), guide="none") +
    scale_x_continuous(breaks=seq.int(1, npc, 1)) +
    labs(title = "Scree plot", subtitle=paste0(nvars, " features, ", npc, " PCs"), 
         x = "Principal Component", y='% Variance') +
    theme(panel.grid=element_blank())
  
  
  # Top influencers. 
  
  infl <- as.data.table(pca$rotation)
  highinfl <- lapply(infl, function(x){ order(abs(x), decreasing=TRUE) })
  infl[, rowID := rownames(pca$rotation)]
  highinfl <- lapply(highinfl, function(x){ head(unique(infl$rowID[x]), ntop) })
  highinfl <- highinfl[1:topnpc]
  setkey(infl, rowID)
  
  # Influence plots
  # lapply(1:npc, function(i) {
  #   ggplot(infl[highinfl[[i]], 
  #             c('rowID', names(loads)[i]), with=FALSE], 
  #        aes_string(x="rowID", xend="rowID", y=0, yend=names(loads)[i])) +
  #   geom_segment() +
  #   geom_point(aes_string(y=names(loads)[i])) +
  #   coord_flip() +
  #   labs(y=NULL, x=NULL)
  # })
  
  # Coefficient names.
  
  ig <- names(covars)
  ig <- ig[! ig %in% c('sample', 'Sample', 'name', 'Name', 'sizeFactor')]
  
  # Means and Variances of the selected genes
  genevar[, selected := name %in% colnames(subcmat)]
  pvar1 <- ggMarginal(
    ggplot(genevar, aes(x=Mean, y=StDev, label=name, colour=selected)) +
      geom_point(shape=16, size=0.8, alpha=0.3) +
      scale_x_log10() +
      scale_y_log10() +
      scale_colour_manual(values=c("black", "red")) +
      annotation_logticks(base=10, sides='lb') +
      labs(x="Mean", y="Standard Deviation", title=paste("Features:", nrow(genevar), "total"), 
           subtitle=paste(nrow(genevar[StDev>0,]), "variable,", length(topgenes), "used")) +
      theme(panel.grid=element_blank()),
    type = "histogram")
  
  # Plot higher PCs in 2D pairs. Highlight one variable at a time.
  L <- lapply(seq(1, topnpc, 2), function(pcx){
    pcy <- pcx + 1
    if (pcy > npc) pcy <- 1
    pca_plotter(pc, infl, highinfl, ig, pcaimp, pcx, pcy)
  } )
  
  L1 <- lapply(L, function(x) { x[[1]] })
  L2 <- lapply(L, function(x) { x[[2]] })
  
  out <- list(pca=pca,
              nvars=nvars,
              nPC=npc,
              pimp=pimp, pvar1=pvar1,
              highloads=highinfl,
              pc2d=L1, pcv2d=L2
  )
  
  if (!is.null(rds)) {
    saveRDS(out, file = rds)
  }
  
  return(out)
}



### GGpairs for continuous variables within a dataframe.
########################################################

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




## GO enrichment
################

do_go <- function(de, db, ntop, onto){
  if (any(de$meetthresh)) {
    geneList <- factor(as.integer(de$meetthresh & de$istop))
    names(geneList) <- de$name

    GOdata <- new("topGOdata",
                  description = paste(onto, db),
                  ontology = onto,
                  allGenes = geneList,
                  nodeSize = 10,
                  annot = annFUN.org,
                  mapping = db,
                  ID = "ensembl")
    resFish <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    allRes <- GenTable(GOdata, classicFisher = resFish, ranksOf = "classicFisher", topNodes = ntop)

    return(allRes[allRes$classicFisher < 0.05,])
  }
}

# for (res  in resall){
#   # GO term enrichment
#   print(paste(res[["contrast"]], "-- GO Biological Process"))
#   go1 <- do_go(res[["shrunkLFC"]], 'org.Mm.eg', ntop, 'BP')
#   print(go1)
#   print(paste(res[["contrast"]], "-- GO Molecular Function"))
#   go2 <- do_go(res[["shrunkLFC"]], 'org.Mm.eg', ntop, 'MF')
#   print(go2)
#   print(paste(res[["contrast"]], "-- GO Cellular Component"))
#   go3 <- do_go(res[["shrunkLFC"]], 'org.Mm.eg', ntop, 'CC')
#   print(go3)
# }
