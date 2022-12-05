


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
