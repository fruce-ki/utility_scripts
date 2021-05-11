


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
    geom_bar(stat='identity') +
    scale_fill_identity() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          legend.position = 'none')
}

# colourblind palette
showpalette( c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") )


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
my_pairwise_internal_corels <- function(mat = rpm, samples = samples_tidy$name, method = "pearson", prefix=params$outpref, txs=4) {
  # Correlations
  cormat <- cor(mat, method=method)
  
  # Cluster
  hcfit <- hclust(dist(scale(cormat, center=TRUE)))
  rn <- rownames(cormat)
  cormat <- cormat[samples, samples]                    # Supplied order
  cormat2 <- cormat                                     # Duplicate on which to delete the diagonal with original order.
  cormat3 <- cormat[rn[hcfit$order], rn[hcfit$order]]   # Duplicate in clustered order.
  cormat4 <- cormat3                                    # duplicate on which to delete the diagonal with in clustered order.
  
  # Delete diagonal half for the numeric labels.
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
  
  rn2 <- rownames(cormat2)
  cormat2 <- as.data.table(cormat2)
  cormat2[, observation1 := factor(rn2, ordered=TRUE, levels=rn2)]
  cormat2 <- melt(cormat2, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat2[, observation2 := factor(observation2, ordered=TRUE, levels=rn2)]
  cormat2 <- cormat2[!is.na(Correlation)]
  
  rn3 <- rownames(cormat3)
  cormat3 <- as.data.table(cormat3)
  cormat3[, observation1 := factor(rn3, ordered=TRUE, levels=rn3)]
  cormat3 <- melt(cormat3, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat3[, observation2 := factor(observation2, ordered=TRUE, levels=rn3)]
  
  rn4 <- rownames(cormat4)
  cormat4 <- as.data.table(cormat4)
  cormat4[, observation1 := factor(rn3, ordered=TRUE, levels=rn4)]
  cormat4 <- melt(cormat4, id.vars = "observation1", value.name = "Correlation", variable.name = "observation2")
  cormat4[, observation2 := factor(observation2, ordered=TRUE, levels=rn4)]
  cormat4 <- cormat4[!is.na(Correlation)]
  
  # Text colour switch for the dynamic range
  m <- min(cormat4$Correlation, na.rm=TRUE)
  M <- max(cormat4$Correlation, na.rm=TRUE)
  colourswitch <- c( m + 0.49 * (M-m),  m + 0.51 * (M-m) )
  
  
  # Square. Custom order. No values. Full range.
  p1 <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Custom order. No values. Dynamic range.
  p1a <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Triangle. Custom order. With values. Full range.
  p2 <- ggplot(cormat2, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Triangle. Custom order. With values. Dynamic range.
  p2a <- ggplot(cormat2, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
    scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Custom order. With values. Full range.
  p12 <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat2, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Custom order. With values. Dyhamic range.
  p12a <- ggplot(cormat, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat2, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
    scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  
  # Square. Clustered order. No values. Full range.
  p3 <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Clustered order. No values. Dyhamic range.
  p3a <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Triangle. Clustered order. with values. Full range.
  p4 <- ggplot(cormat4, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Triangle. Clustered order. with values. Dyhamic range.
  p4a <- ggplot(cormat4, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "transparent" ) +
    scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Clustered order. With values. Full range.
  p34 <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat4, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=Correlation >= -0.60 & Correlation <= 0.60 ), size=rel(txs)) +
    scale_fill_gradientn(limits=c(-1, 1), colors=c("lightskyblue", "dodgerblue3", "darkblue", "black", "darkred", "red", "gold"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c('FALSE'="black", 'TRUE'="white"), na.value="forestgreen", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  # Square. Clustered order. With values. Dyhamic range.
  p34a <- ggplot(cormat3, aes(x=observation1, y=observation2)) +
    geom_tile(aes(fill=Correlation)) +
    geom_text(data=cormat4, aes(label=sub('0.', '.', as.character(round(Correlation, 2))), colour=(Correlation <= colourswitch[2]) ), size=rel(txs)) +
    scale_fill_gradientn(colors=c("black", "red", "gold", "white"), na.value = "forestgreen" ) +
    scale_colour_manual(values=c("black", "white"), na.value="transparent", guide="none") +
    scale_x_discrete(position = "top") +
    labs(x='', y='', title=paste(paste(toupper(substr(method, 1, 1)), tolower(substr(method, 2, nchar(method))), "'s", sep=""), "correlation - Clustered")) +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
          panel.grid = element_blank() )
  
  fwrite(dcast(cormat2, observation1 ~ observation2, value.var = "Correlation"),
         file=paste0(prefix, '_cor.txt'),
         sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  return( list(sfrnc=p1, tfrnc=p2, frnc=p12,
               sdrnc=p1a, tdrnc=p2a, drnc=p12a,
               sfrc=p3, tfrc=p4, frc=p34,
               sdrc=p3a, tdrc=p4a, drc=p34a) )
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


### PCA and plots
#################
# covars is DE style colData
do_pca <- function(countsmat, covars, center = TRUE, scale = TRUE, topvars=params$topvars, loadthresh = params$minL, prefix = paste0(params$prefix, '_pc')) {
  # countsmat = counts
  
  # Gene variances
  genevar <- data.table(name = rownames(countsmat),
                        Mean = rowMeans(countsmat),
                        Var = rowVars(countsmat) )
  
  # Use only the most variable genes, and ensure they are variable.
  setorder(genevar, -Var)
  topgenes = genevar[1:min(topvars, nrow(genevar)), ] [Var>0, name]
  countsmat <- countsmat[topgenes,]
  
  # Rotate counts so genes are variables and samples are observations
  countsmat <- t(countsmat)
  nvars <- dim(countsmat)[2]
  message(paste("Number of most variable features considered:", nvars))
  
  # Add covariates info
  # counts <- cbind(covars, as.data.table(countsmat))
  
  pca <- prcomp(countsmat, center = center, scale. = scale)
  srn <- sqrt(nrow(countsmat) - 1)
  
  # pc <- sweep(pca$x, 2, 1 / (pca$sdev * srn), FUN = '*')   # wrong?
  pc <- pca$x
  
  pc <- cbind(pc, data.frame(sample = rownames(pc)))
  #covars <- cbind(covars, data.frame(sample = rownames(covars)))
  pc <- as.data.frame(merge(pc, covars, by="sample", all = TRUE))
  npc <- sum(pca$sdev > 1)
  
  # Screeplot.
  
  pcaimp <- as.data.table(cbind(as.data.frame(t(summary(pca)$importance)),
                                data.frame(PC = 1:length(colnames(pca$x)))))
  pcaimp <- melt(pcaimp, variable.name = "type", value.name = "Proportion", id.vars = c("PC", "Standard deviation"))
  
  pimp <- ggplot(pcaimp, aes(x=PC, y=Proportion)) +
    facet_grid(type ~ ., scales="free_y") +
    geom_bar(aes(fill=type), stat="identity", width=0.5) +
    geom_text(aes(x=PC-0.3, y=0, label=paste0(round(Proportion*100, 1),"%"), colour=rev(type)), angle=90, hjust=0, vjust=0, size=rel(3)) +
    geom_text_repel(aes(label=paste0(round(Proportion*100, 1),"%"), colour=type), ylim=c(0, max(pcaimp[type=='Cumulative Proportion', Proportion])), angle=90, vjust=0, direction= 'y', size=rel(3)) +
    geom_line(aes(colour=type)) +
    scale_fill_manual(values = c("grey50", "transparent")) +
    scale_colour_manual(values = c("transparent", "grey25")) +
    labs(title = "Scree plot", subtitle=paste0(nvars, " features, ", npc, " PCs"), x = "Principal Component") +
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
  
  pvar1 <- ggMarginal(ggplot(genevar, aes(x=Mean, y=StDev, label=rowID, colour=istop)) +
                        geom_point(shape=16, size=0.8, alpha=0.5) +
                        scale_x_log10() +
                        scale_colour_manual(values=c("black", "violetred")) +
                        scale_y_log10() +
                        labs(x="Mean", y="Standard Deviation", title=paste0("All features (", nvars, ")")) +
                        theme(legend.position = "none"),
                      type = "histogram")
  
  pvar2 <- ggplot(genevar[(istop),], aes(x=Mean, y=StDev, label=rowID)) +
    geom_point(shape=16, size=0.8, alpha=0.5, colour="purple") +
    labs(x="Mean", y="Standard Deviation",
         title=paste0("High-loading features (", length(sel), ")"), subtitle = paste0("( >", loadthresh, ") in any of the first 3 PCs")) +
    scale_x_log10() +
    scale_y_log10()
  
  # Correlation of the samples for the selected features only.
  
  pcor <- my_pairwise_internal_corels(t(countsmat[, sel]), samples=covars$sample, prefix=prefix)
  
  # Plot first 3 PCs in 2D pairs.
  # Highlight one variable at a time.
  pc12 <- lapply(ig, function(varname) {
    return(
      ggplot(pc, aes_string(x="PC1", y="PC2", label="sample", colour=varname)) +
        geom_point(alpha=0.8, size=rel(1.5)) +
        coord_fixed() +
        labs(x=paste0("PC1 (", round(pcaimp[1,4]*100, 1), "%)"), y=paste0("PC2 (", round(pcaimp[2,4]*100, 1), "%)"))
    )})
  
  pc13 <- lapply(ig, function(varname) {
    return(
      ggplot(pc, aes_string(x="PC1", y="PC3", label="sample", colour=varname)) +
        geom_point(alpha=0.8, size=rel(1.5)) +
        coord_fixed() +
        labs(x=paste0("PC1 (", round(pcaimp[1,4]*100, 1), "%)"), y=paste0("PC3 (", round(pcaimp[3,4]*100, 1), "%)"))
    )})
  
  pc32 <- lapply(ig, function(varname) {
    return(
      ggplot(pc, aes_string(x="PC3", y="PC2", label="sample", colour=varname)) +
        geom_point(alpha=0.8, size=rel(1.5)) +
        coord_fixed() +
        labs(x=paste0("PC3 (", round(pcaimp[3,4]*100, 1), "%)"), y=paste0("PC2 (", round(pcaimp[2,4]*100, 1), "%)"))
    )})
  
  return(list(pca=pca,
              nvars=nvars,
              nPC=npc,
              pimp=pimp,
              load1=sel1, load2=sel2, load3=sel3,
              pc_1_2=pc12, pc_1_3=pc13, pc_3_2=pc32,
              pvar1=pvar1, pvar2=pvar2,
              pcor1=pcor[['frnc']], pcor2=pcor[['drc']]
  ))
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
