### PCA
########

library(data.table)
library(matrixStats)
library(ggplot2)
library(ggExtra)
library(plotly)
library(htmltools)
# Calculate PCA and plots
do_pca <- function(counts, covars, center=TRUE, scale=TRUE, loadthresh = 0.75){
  # Convert table to matrix
  countsmat <- as.matrix(counts[, 2:length(counts)])
  rownames(countsmat) <- counts$Geneid
  
  # Gene variances
  genevar <- data.table(name = rownames(countsmat),
                        Mean = rowMeans(countsmat),
                        StDev = rowVars(countsmat) )
  
  # Rotate counts so genes are variables and samples are observations
  countsmat <- t(countsmat)
  
  # Add covariates info
  counts <- cbind(covars, as.data.table(countsmat))
  rownames(counts) <- rownames(countsmat)
  
  # Calculate principal components (after getting rid of zero-variance genes that can't be standardized).
  
  nonconstant <- which(apply(countsmat, 2, var)!=0)
  countsmat <- countsmat[, nonconstant]
  nvars <- dim(countsmat)[2]
  
  pca <- prcomp(countsmat, center = center, scale. = scale)
  srn <- sqrt(nrow(countsmat) - 1)
  pc <- sweep(pca$x, 2, 1 / (pca$sdev * srn), FUN = '*')
  pc <- cbind(pc, data.frame(name = rownames(pc)))
  pc <- merge(pc, covars, by.x = "name", by.y="sample", all = TRUE)
  npc <- sum(pca$sdev > 1)
  
  # Screeplot.
  
  pcaimp <- as.data.table(cbind(as.data.frame(t(summary(pca)$importance)), 
                                data.frame(PC = 1:length(colnames(pca$x)))))
  pcaimp <- melt(pcaimp, variable.name = "type", value.name = "Proportion", id.vars = c("PC", "Standard deviation"))
  
  pcmax=15
  pimp <- ggplot(pcaimp[(PC<pcmax),], aes(x=PC, y=Proportion)) +
    facet_grid(type ~ ., scales="free_y") +
    geom_bar(aes(fill=type), stat="identity", width=0.6) +
    geom_line(aes(colour=type)) +
    scale_fill_manual(values = c("grey50", "transparent")) +
    scale_colour_manual(values = c("transparent", "grey25")) +
    scale_x_continuous(breaks = seq(1, pcmax, 2)) +
    labs(title = "PCA Screeplot", subtitle=paste0(nvars, " features, ", npc, " PCs"), x = "Principal Component") +
    theme(axis.title.y = element_blank(),
          legend.position = "none")
  
  # Top loadings for the first 3 PCs.
  
  loads <- as.data.table(t(t(pca$rotation) * pca$sdev))
  loads[, rowID := rownames(pca$rotation)]
  loads <- melt(loads, value.name = "Loading", id.vars="rowID", variable.name="PC")
  loads[, absload := abs(Loading)]
  lhigh <- loads$absload > loadthresh     # keep only high loadings
  setorder(loads, PC, -absload)
  sel1 <- loads[lhigh > loadthresh & PC=="PC1", rowID]
  sel2 <- loads[lhigh > loadthresh & PC=="PC2", rowID]
  sel3 <- loads[lhigh > loadthresh & PC=="PC3", rowID]
  sel <- unique(loads[lhigh, rowID])
  
  # Coefficient names.
  
  ig <- names(covars)
  ig <- ig[! ig %in% c('sample', 'Sample', 'name', 'Name', 'sizeFactor')]
  
  # Means and Variances of the top loads in the first 3 PCs.
  
  genevar = data.table(rowID = colnames(countsmat),
                       Mean = colMeans(countsmat), 
                       StDev = colSds(countsmat) )
  genevar[, istop:= rowID %in% sel]
  
  pvar1 <- ggMarginal(ggplot(genevar) +
                        geom_point(aes(x=Mean, y=StDev, label=rowID, colour=istop), shape=16, size=0.8, alpha=0.5) +
                        scale_x_log10() +
                        scale_colour_manual(values=c("black", "violetred")) +
                        scale_y_log10() +
                        labs(x="Mean", y="Standard Deviation", title="All features") +
                        theme(legend.position = "none"),
                      type = "histogram")
  
  pvar2 <- ggplot(genevar[(istop),]) +
    geom_point(aes(x=Mean, y=StDev, label=rowid), shape=16, size=0.8, alpha=0.5, colour="purple") +
    labs(x="Mean", y="Standard Deviation", 
         title="Top-loading features", subtitle = "in the first 3 PCs") +
    scale_x_log10() +
    scale_y_log10()
  
  # Correlation of the samples for the selected features only.
  
  genec <- cor(t(countsmat[, sel]))
  hcfit <- hclust(dist(t(genec)))
  rn <- rownames(genec)
  genec <- as.data.table(genec)
  genec[, observation1 := rn]
  genec <- melt(genec, id.vars = "observation1", value.name = "correlation", variable.name = "observation2")
  genec[, observation1 := factor(observation1, ordered=TRUE, 
                                 levels=rn[hcfit$order]) ]
  genec[, observation2 := factor(observation2, ordered=TRUE, 
                                 levels=rn[hcfit$order]) ]
  
  pcor1 <- ggplot(genec, aes(x=observation1, y=observation2, fill=correlation, label=correlation)) +
    geom_tile() +
    # geom_text(colour="white") +
    # scale_x_discrete(position="top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("steelblue2", "blue", "darkblue", "black", "darkred", "red", "orange", "yellow") ) +
    labs(x='', y='', title="Correlation of samples", subtitle="(top loading in first 3PCs)") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  
  # # Correlation of the selected features.
  # 
  # genec <- cor(countsmat[, sel])
  # hcfit <- hclust(dist(t(genec)))
  # rn <- rownames(genec)
  # genec <- as.data.table(genec)
  # genec[, observation1 := rn]
  # genec <- melt(genec, id.vars = "observation1", value.name = "correlation", variable.name = "observation2")
  # genec[, observation1 := factor(observation1, ordered=TRUE, 
  #                                levels=rn[hcfit$order]) ]
  # genec[, observation2 := factor(observation2, ordered=TRUE, 
  #                                levels=rn[hcfit$order]) ]
  # 
  # pcor2 <- ggplot(genec, aes(x=observation1, y=observation2, fill=correlation, label=correlation)) +
  #   geom_tile() +
  #   # geom_text(colour="white") +
  #   # scale_x_discrete(position="top") +
  #   scale_fill_gradientn(limits=c(-1, 1), colors=c("steelblue2", "blue", "darkblue", "black", "darkred", "red", "orange") ) +
  #   labs(x='', y='', title="Correlation of top-loading features", subtitle="(in first 3 PCs)") +
  #   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  # if (sqrt(dim(genec)[1]) >50) {
  #   pcor2 <- pcor2 + theme(axis.text.x = element_blank(), 
  #                          axis.text.y = element_blank(),
  #                          axis.ticks = element_blank() )
  # }
  
  # Plot first 3 PCs in 2D pairs.
  # Highlight one variable at a time.
  pc12 <- lapply(ig, function(varname) { 
    return( 
      ggplot(pc, aes_string(x="PC1", y="PC2", label="name", colour=varname)) +
        geom_point(shape=16) +
        coord_fixed()
    )})
  # do.call(grid.arrange, c(l1, ncol=2))
  
  pc13 <- lapply(ig, function(varname) { 
    return( 
      ggplot(pc, aes_string(x="PC1", y="PC3", label="name", colour=varname)) +
        geom_point(shape=16) +
        coord_fixed()
    )})
  # do.call(grid.arrange, c(l2, ncol=2))
  
  pc23 <- lapply(ig, function(varname) { 
    return( 
      ggplot(pc, aes_string(x="PC2", y="PC3", label="name", colour=varname)) +
        geom_point(shape=16) +
        coord_fixed()
    )})
  # do.call(grid.arrange, c(l3, ncol=2))
  
  return(list(pca=pca, 
              nvars=nvars, 
              nPC=npc, 
              pimp=pimp, 
              load1=sel1, load2=sel2, load3=sel3,
              pc_1_2=pc12, pc_1_3=pc13, pc_2_3=pc23,
              pvar1=pvar1, pvar2=pvar2, 
              pcor1=pcor1 #, pcor2=pcor2
  ))
}

### Correlations within the dataframe.
######################################

# library(ggplot2)
# library(data.table)
# colanes and rownames yes, non-numeric columns no.
my_pairwise_internal_coords <- function(mat, method="pearson", triangular=FALSE) {
  # Correlations
  cormat <- cor(mat, method=method)
  
  # Cluster
  hcfit <- hclust(dist(scale(cormat, center=TRUE)))
  cormat <- cormat[hcfit$order, hcfit$order]
  rn <- rownames(cormat)
  cormat2 <- cormat # Duplicate,
  
  
  # Delete diagonal half for the numeric labels. But keep the heatmap square, because pretty.
  for (r in 1:nrow(cormat)) {
    for (c in 1:ncol(cormat)) {
      if (c < r) {
        cormat[r, c] <- NA_real_
      }
    }
  }
  # Unless I want it to be all triangular
  if (triangular) {
    for (r in 1:nrow(cormat2)) {
      for (c in 1:ncol(cormat2)) {
        if (c < r) {
          cormat2[r, c] <- NA_real_
        }
      }
    }
  }
  
  cormat <- as.data.table(t(cormat))
  cormat[, observation1 := factor(rn, ordered=TRUE, levels=rn)]
  cormat <- melt(cormat, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat[, observation2 := factor(observation2, ordered=TRUE, levels=rn)]
  cormat2 <- as.data.table(t(cormat2))
  cormat2[, observation1 := factor(rn, ordered=TRUE, levels=rn)]
  cormat2 <- melt(cormat2, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat2[, observation2 := factor(observation2, ordered=TRUE, levels=rn)]
  
  p <- ggplot(cormat2, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat, aes(label=round(Correlation, 2), colour=Correlation)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue", "blue", "darkblue", "black", "darkred", "red", "orange", "yellow"), na.value = "transparent" ) +
    scale_colour_gradientn(limits=c(-1, 1), colors=rep(c("black", "white", "white", "black"), each=3), guide="none") +
    labs(x='', y='') +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  if (!triangular) {
    p <- p + theme(panel.grid.major = element_blank())
  }
  return(p)
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
