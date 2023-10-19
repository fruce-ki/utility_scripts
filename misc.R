
# RMarkdown header template
###########################

# title: "Report"
# author: "Kimon Froussios"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
# output:
#   html_document:
#     code_folding: hide
#     toc:
#       toc_float: true
#       toc_depth: 4
#       collapsed: false
# params:




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
  intensity_colours <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  showpalette(intensity_colours)
  
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
{
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
  
  ### Pairwise correlations of all columns of a matrix
  ####################################
  my_pairwise_internal_correls <- function(mat, samples = NULL, method = "pearson", minMean=0, minSingle=0) {
    # mat <- params$tpm[topgenes,]; samples <- params$covars$Sample
    # samples = NULL
    if (is.null(samples))
      samples <- colnames(mat)
    
    # Filter
    if (minMean != 0 | minSingle != 0) {
      mat <- mat[rowSums(mat >= minSingle) >= 1 | rowMeans(mat) >= minMean, ]
      cat("Minimum mean expression across all samples:", minMean, "\n")
      cat("Minimum expression in any single sample, if the minimum mean is not satsfied:", minSingle, "\n")
    } else {
      mat <- mat[rowSums(mat) > 0, ]
    }
    cat("Number of eligible features:", nrow(mat), "\n")
    
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
    
    return(list(unord = cormat, unordtri = cormat2, clust = cormat3, clusttri = cormat4, meth = method, dendroR = pd, dendroC = NULL))
  }
  
  
  ### Pairwise correlations of all columns between two matrices
  ####################################
  my_pairwise_external_correls <- function(a, b, method = "pearson", minMean=0, minSingle=0, logged = FALSE) {
    # a = BSP; b = NLT; method = "pearson"; minMean = 0; minSingle = 0; logged = TRUE
    if ( ! all(rownames(a) == rownames(b)) )
      stop("The rownames don't match up between the datasets.")
    
    # Combine into a single matrix, then use the internal_correls function, then delete rows and columns that correspond to internal correlations of each subset.
    # Computationally more intensive than it needs to be, but simpler to implement given the pre-existing code, and easier to keep consistent.
    # mat <- cbind(a, b)
    # mycors <- my_pairwise_internal_correls(mat, method=method, minMean=minMean, minSingle=minSingle)
    # 
    # mycors[1:4] <- lapply(mycors[1:4], function(z){
    #   # z = mycors[[2]]
    #   z[rownames(z) %in% colnames(a), colnames(z) %in% colnames(b)]
    # })
    # 
    # return(mycors)
    
    # The above method skews the clustering of samples. It does not allow the axes to cluster independently.
    # Re-implement fully.
    
    # Filter. Both datasets must come out with the same surviving rows, in order to remain comparable.
    if (minMean != 0 | minSingle != 0) {
      sel <- (rowSums(a >= minSingle) >= 1 | rowMeans(a) >= minMean) |
             (rowSums(b >= minSingle) >= 1 | rowMeans(b) >= minMean)
      cat("Minimum mean expression across all samples in each dataset:", minMean, "\n")
      cat("Minimum expression in any single sample in each dataset, if the minimum mean is not satisfied:", minSingle, "\n")
    } else {
      sel <- rowSums(a) > 0 | 
             rowSums(b) > 0
    }
    a <- a[sel, ]
    b <- b[sel, ]
    cat("Number of eligible features:", nrow(a), "\n")
    
    if (logged) {
      # Handle zeros
      pad <- min(c(a[a>0], b[b>0])) / 10 
      
      a <- log10(a)
      b <- log10(b)
    }
    
    # Correlations
    cormat <- cor(a, b, method=method)
    
    # Cluster rows
    hcfitr <- hclust(dist(scale(cormat, center=TRUE)))
    rn <- rownames(cormat)
    # Make dendrogram. https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
    dendr <- as.dendrogram(hcfitr)
    dendr_data <- dendro_data(dendr)
    # Setup the data, so that the axes are exchanged, instead of using coord_flip()
    segmentr_data <- with(
        segment(dendr_data), 
        data.frame(x = y, y = x, xend = yend, yend = xend))
    # Use the dendrogram label data to position the sample labels
    sample_posr_table <- with(
        dendr_data$labels, 
        data.frame(y_center = x, Sample = as.character(label), height = 1))
    # Limits for the vertical axis
    sample_axisr_limits <- with(
        sample_posr_table, 
        c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
    ) + 0.1 * c(-1, 1) # extra spacing: 0.1
    # Dendrogram plot
    pdr <- ggplot(segmentr_data) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      scale_x_reverse(expand = c(0, 0.5),
                      position = "top") + 
      scale_y_continuous(position = "right",
                         breaks = sample_posr_table$y_center, 
                         labels = sample_posr_table$Sample, 
                         limits = sample_axisr_limits, 
                         expand = c(0, 0)) + 
      labs(x = NULL, y = NULL) +
      theme_minimal() + 
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  
    # Cluster columns
    hcfitc <- hclust(dist(scale(t(cormat), center=TRUE)))
    cn <- colnames(cormat)
    # Make dendrogram.
    dendc <- as.dendrogram(hcfitc)
    dendc_data <- dendro_data(dendc)
    segmentc_data <- with(
        segment(dendc_data), 
        data.frame(x = x, y = y, xend = xend, yend = yend))
    sample_posc_table <- with(
        dendc_data$labels, 
        data.frame(x_center = x, Sample = as.character(label), height = 1))
    sample_axisc_limits <- with(
        sample_posc_table, 
        c(min(x_center - 0.5 * height), max(x_center + 0.5 * height))
    ) + 0.1 * c(-1, 1) # extra spacing: 0.1
    pdc <- ggplot(segmentc_data) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      scale_y_continuous(expand = c(0, 0.5),
                         position = "left") + 
      scale_x_continuous(position = "bottom",
                         breaks = sample_posc_table$x_center, 
                         labels = sample_posc_table$Sample, 
                         limits = sample_axisc_limits, 
                         expand = c(0, 0)) + 
      labs(x = NULL, y = NULL) +
      theme_minimal() + 
      theme(panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  
    # Create duplicates for different plot styles
    cormat2 <- cormat                                     # Duplicate in which to delete below the diagonal.
    cormat3 <- cormat[rn[hcfitr$order], cn[hcfitc$order]]   # Duplicate in clustered order.
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
    
    return(list(unord = cormat, unordtri = cormat2, clust = cormat3, clusttri = cormat4, meth = method, dendroR = pdr, dendroC = pdc))
  }
  
  
  ### Plot the results of the internal or external pairwise correlations
  ####################################
  plot_my_correlations <- function(matlist, rds=NULL, txs=3) {
    # matlist <- my_pairwise_internal_correls(rawBSP); rds <- NULL; txs <- 3
    # Extract. Not necessary, but quicker than refactoring the code after separating it from my_pairwise_internal_correls
    cormat <- matlist[["unord"]]
    cormat2 <- matlist[["unordtri"]]
    cormat3 <- matlist[["clust"]]
    cormat4 <- matlist[["clusttri"]]
    method <- matlist[["meth"]]
    denr <- matlist[["dendroR"]]
    denc <- matlist[["dendroC"]]
  
    # Restructure for plotting.
    rn <- rownames(cormat)
    cn <- colnames(cormat)
    cormat <- as.data.table(cormat)
    cormat[, observation1 := factor(rn, ordered=TRUE, levels=rn)]
    cormat <- melt(cormat, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
    cormat[, observation2 := factor(observation2, ordered=TRUE, levels=cn)]
    # cormat <- merge(cormat, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
    
    rn2 <- rownames(cormat2)
    cn2 <- colnames(cormat2)
    cormat2 <- as.data.table(cormat2)
    cormat2[, observation1 := factor(rn2, ordered=TRUE, levels=rn2)]
    cormat2 <- melt(cormat2, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
    cormat2[, observation2 := factor(observation2, ordered=TRUE, levels=cn2)]
    cormat2 <- cormat2[!is.na(Correlation)]
    # cormat2 <- merge(cormat2, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
    
    rn3 <- rownames(cormat3)
    cn3 <- colnames(cormat3)
    cormat3 <- as.data.table(cormat3)
    cormat3[, observation1 := factor(rn3, ordered=TRUE, levels=rn3)]
    cormat3 <- melt(cormat3, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
    cormat3[, observation2 := factor(observation2, ordered=TRUE, levels=cn3)]
    # cormat3 <- merge(cormat3, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
    
    rn4 <- rownames(cormat4)
    cn4 <- colnames(cormat4)
    cormat4 <- as.data.table(cormat4)
    cormat4[, observation1 := factor(rn3, ordered=TRUE, levels=rn4)]
    cormat4 <- melt(cormat4, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
    cormat4[, observation2 := factor(observation2, ordered=TRUE, levels=cn4)]
    cormat4 <- cormat4[!is.na(Correlation)]
    # cormat4 <- merge(cormat4, sample_pos_table, by.x="observation2", by.y="Sample", all.x=TRUE)
    
    # Text colour switch for the dynamic range
    m <- min(cormat4$Correlation, na.rm=TRUE)
    M <- max(cormat4$Correlation, na.rm=TRUE)
    colourswitch <- c( m + 0.49 * (M-m),  m + 0.51 * (M-m) )
    
    
    # Square. Custom order. No values. Full range.
    pfr <- ggplot(cormat, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
      scale_x_discrete(position = "top") +
      labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
  
    # Square. Custom order. No values. Dynamic range.
    pdr <- ggplot(cormat, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
      scale_x_discrete(position = "top") +
      labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
  
    # # Triangle. Custom order. With values. Full range.
    # pfrt <- ggplot(cormat2, aes(y=observation1, x=observation2)) +
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
    # pdrt <- ggplot(cormat2, aes(y=observation1, x=observation2)) +
    #   geom_tile(aes(fill=Correlation)) +
    #   geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    #   scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
    #   scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    #   scale_x_discrete(position = "top") +
    #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    #   theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
    #         panel.grid = element_blank() )
  
    # Square. Custom order. With values triangle. Full range.
    pfrv <- ggplot(cormat, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      geom_text(data=cormat2, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
      scale_x_discrete(position = "top") +
      scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
      scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
      labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
  
    # Square. Custom order. With values triangle. Dynamic range.
    pdrv <- ggplot(cormat, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      geom_text(data=cormat2, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
      scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
      scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
      scale_x_discrete(position = "top") +
      labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
  
    
    # Square. Clustered order. No values. Full range.
    pfrc <- ggplot(cormat3, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
      scale_x_discrete(position = "top") +
      labs(x='', y='', caption=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
  
    # Square. Clustered order. No values. Dynamic range.
    pdrc <- ggplot(cormat3, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
      scale_x_discrete(position = "top") +
      labs(x='', y='', caption=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
  
    # # Triangle. Clustered order. With values. Full range.
    # pfrtc <- ggplot(cormat4, aes(y=observation1, x=observation2)) +
    #   geom_tile(aes(fill=Correlation)) +
    #   geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    #   scale_x_discrete(position = "top") +
    #   scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    #   scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
    #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    #   theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
    #         panel.grid = element_blank() )
    # 
    # # Triangle. Clustered order. With values. Dynamic range.
    # pdrtc <- ggplot(cormat4, aes(y=observation1, x=observation2)) +
    #   geom_tile(aes(fill=Correlation)) +
    #   geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    #   scale_x_discrete(position = "top") +
    #   scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
    #   scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    #   labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    #   theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
    #         panel.grid = element_blank() )
  
    # Square. Clustered order. With values triangle. Full range.
    pfrvc <- ggplot(cormat3, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      geom_text(data=cormat4, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
      scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
      scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
      scale_x_discrete(position = "top") +
      labs(x='', y='', caption=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
  
    # Square. Clustered order. With values triangle. Dynamic range.
    pdrvc <- ggplot(cormat3, aes(y=observation1, x=observation2)) +
      geom_tile(aes(fill=Correlation)) +
      geom_text(data=cormat4, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
      scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
      scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
      scale_x_discrete(position = "top") +
      labs(x='', y='', caption=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
      theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
            panel.grid = element_blank() )
    
    if(is.null(denc)){
      out <- list(corr=dcast(cormat2, observation1 ~ observation2, value.var = "Correlation"),
                  pfr = pfr, pdr = pdr,
                  # pfrt = pfrt, pdrt = pdrt,
                  pfrv = pfrv, pdrv = pdrv,
                  pfrc = denr + pfrc + plot_layout(ncol=2, nrow=1, widths=c(1,4)),
                  pdrc = denr + pdrc + plot_layout(ncol=2, nrow=1, widths=c(1,4)),
                  # pfrtc = pfrtc, pdrtc = pdrtc,
                  pfrvc = denr + pfrvc + plot_layout(ncol=2, nrow=1, widths=c(1,4)),
                  pdrvc = denr + pdrvc + plot_layout(ncol=2, nrow=1, widths=c(1,4))
                  )
    } else {
      out <- list(corr=dcast(cormat2, observation1 ~ observation2, value.var = "Correlation"),
                  pfr = pfr, pdr = pdr,
                  # pfrt = pfrt, pdrt = pdrt,
                  pfrv = pfrv, pdrv = pdrv,
                  pfrc = plot_spacer() + denc + denr + pfrc + plot_layout(ncol=2, nrow=2, widths=c(1,4), heights = c(1,7)), 
                  pdrc = plot_spacer() + denc + denr + pdrc + plot_layout(ncol=2, nrow=2, widths=c(1,4), heights = c(1,7)),
                  # pfrtc = pfrtc, pdrtc = pdrtc,
                  pfrvc = plot_spacer() + denc + denr + pfrvc + plot_layout(ncol=2, nrow=2, widths=c(1,4), heights = c(1,7)),
                  pdrvc = plot_spacer() + denc + denr + pdrvc + plot_layout(ncol=2, nrow=2, widths=c(1,4), heights = c(1,7))
                  )
    }
    
    if (!is.null(rds) && length(rds) > 0) {
      saveRDS(out, file = rds)
    }
    
    return(out)
  }

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
  

## GO enrichment
################
{
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
  
}
