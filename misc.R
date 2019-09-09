### Correlations within the dataframe.
######################################
my_pairwise_internal_coords <- function(df) {
  # Correlations
  cm <- cor(df)
  # Order of columns and rows
  cm <- cm[names(df), names(df)]
  # Drop duplicates below the diagonal
  for (x in 1:ncol(cm)){
    for (y in 1:nrow(cm)) {
      if(x < y)
        cm[y, x] <- NA
    }
  }
  # Tidy
  cm <- cm %>%
    as.data.frame() %>%
    rownames_to_column(var="x") %>%
    gather(key="y", value="Correlation", -x) %>%
    mutate(x = factor(x, levels=names(counts[, -c("id")])), y = factor(y, levels=names(counts[, -c("id")])))
  # Plot
  p <- ggplot(cm, aes(x=x, y=y, fill=Correlation, label=round(Correlation, 2))) +
    geom_tile() +
    geom_text(colour="white") +
    scale_x_discrete(position="top") +
    scale_fill_gradient(na.value="transparent") +
    labs(x='', y='') +
    theme(axis.text.x=element_text(angle=45, hjust=0, vjust=0.5),
          panel.background = element_rect(fill="white"))
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
