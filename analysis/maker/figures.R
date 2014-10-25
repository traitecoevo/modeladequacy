fig_cols <- function() {
  c("#F46D43", "#3288BD", "#CDCD00")
}

fig_model_support_ic <- function(fits) {
  type <- guess_ic(fits)
  col <- fig_cols()
  ylab <- paste(toupper(type), "weight")

  ic <- build_ic(fits)

  ## add dummy variable
  dd <- cbind(Subclade=as.character(seq_len(nrow(ic))), ic)
  ord <- dd[order(dd[,"OU"], dd[,"BM"], decreasing = TRUE), "Subclade"]
  ## reorient data frame for plotting
  df <- suppressMessages(melt(dd))
  colnames(df)[2:3] <- c("model", "weight")
  df$Subclade <- factor(df$Subclade, ord)

  ## create geom_bar plot
  .e <- environment()
  p <- ggplot(df, aes(factor(Subclade), weight, fill=model, order=model),
              environment=.e)
  p <- p + geom_bar(stat="identity", position="stack", width=1)
  p <- p + scale_y_continuous(name=ylab)
  p <- p + scale_fill_manual("Model", values=col)
  p <- p + theme_bw()
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text.y=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.border=element_blank(),
                 panel.background=element_blank(),
                 strip.background=element_rect(fill="white"),
                 plot.background=element_blank())
  p <- p + xlab("Dataset")
  ## for now at least, we're going to need to actually run the ggplot
  ## figures.
  print(p)
}

fig_pval_histogram <- function(best) {
  ## prune out irrelevant categories
  keep <- c("trait", pvalue_names_arbutus())
  best <- cbind(rownames(best), best[, c(keep)])

  df <- suppressMessages(melt(best))

  ## set the order and the labels
  df$trait <- rename_traits(df$trait)
  df$variable <- rename_variables(df$variable)

  .e <- environment()

  p <- ggplot(df, aes(x=value), environment = .e)
  p <- p + geom_histogram(binwidth=0.025, alpha=0.8, aes(y=..density.., fill=factor(trait)))
  p <- p + scale_fill_manual(values=fig_cols())
  p <- p + theme_bw()
  p <- p + xlab("p-value")
  p <- p + ylab("Density")
  p <- p + facet_grid(trait~variable, labeller = label_parsed)
  p <- p + theme(strip.background=element_rect(fill="white"),
                 plot.background=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 axis.text=element_text(size=8),
                 axis.ticks.y=element_blank(),
                 axis.text.y=element_blank(),
                 legend.position="none")
  print(p)
}

## Code for plotting the relative support (AIC) for the best model
## (compared to BM) vs. a multivariate measure of model adequacy
fig_modelad_ic <- function(fits) {
  df <- build_table_adequacy_ic(fits)
  df <- prepare_df_for_ggplot(df)

  type <- guess_ic(fits)
  xlab <- sprintf("%s(BM) - %s(OU/EB)", toupper(type), toupper(type))
  col <- fig_cols()
  
  .e <- environment()

  ## need to set options for scientific notation
  options(scipen=1000)
  
  ## the occasional dataset may have NA for Mahalanobis distance
  ## remove this for the plot
  df <- na.omit(df)

  p <- ggplot(df, aes(diff.bm, mv), environment=.e)
  p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.8)
  p <- p + scale_colour_manual("Trait", values=col)
  p <- p + scale_shape_manual("Rank", values=c(15,16,17))
  p <- p + theme_bw()
  p <- p + theme(plot.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border=element_blank(),
                 legend.justification=c(0.02,0.98),
                 legend.position=c(0.02,0.98),
                 legend.key=element_blank(),
                 axis.line = element_line(color = 'black'))
  p <- p + scale_y_log10()
  p <- p + scale_x_log10()
  p <- p + xlab(xlab)
  p <- p + ylab("Mahalanobis distance")
  print(p)
}

## Code to plot a multivariate measure of model adequacy (Mahalanobis
## distance) against clade size
fig_modelad_size <- function(fits_best) {
  df <- prepare_df_for_ggplot(fits_best)
  col <- fig_cols()

  .e <- environment()

  ## the occasional dataset may have NA for Mahalanobis distance
  ## remove this for the plot
  df <- na.omit(df)

  p <- ggplot(df, aes(size, mv), environment=.e)
  p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.8)

  p <- p + scale_colour_manual("Trait", values=col)
  p <- p + scale_shape_manual("Rank", values=c(15,16,17))
  p <- p + theme_bw()
  p <- p + theme(plot.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border=element_blank(),
                 legend.justification=c(0.98,0.02),
                 legend.position=c(0.98,0.02),
                 legend.key=element_blank(),
                 axis.line = element_line(color = 'black'))
  p <- p + scale_y_log10()
  p <- p + scale_x_log10()
  p <- p + xlab("Number of taxa")
  p <- p + ylab("Mahalanobis distance between observed and simulated")
  print(p)
}

## ### Code to plot a multivariate measure of model adequacy (Mahalanobis distance) against clade age
fig_modelad_age <- function(fits_best) {
  df <- prepare_df_for_ggplot(fits_best)
  col <- fig_cols()
  
  .e <- environment()

  ## the occasional dataset may have NA for Mahalanobis distance
  ## remove this for the plot
  df <- na.omit(df)

  p <- ggplot(df, aes(age, mv), environment=.e)
  p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.8)
  p <- p + scale_colour_manual("Trait", values=col)
  p <- p + scale_shape_manual("Rank", values=c(15,16,17))
  p <- p + theme_bw()
  p <- p + theme(plot.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border=element_blank(),
                 legend.justification=c(0.02,0.98),
                 legend.position=c(0.02,0.98),
                 legend.key=element_blank(),
                 axis.line = element_line(color = 'black'))
  p <- p + scale_y_log10()
  p <- p + scale_x_log10()
  p <- p + xlab("Age of crown group (my)")
  p <- p + ylab("Mahalanobis distance")
  print(p)
}

fig_two_clades <- function(examples) {
  me_dat <- examples$Meliaceae$OU
  fa_dat <- examples$Fagaceae$OU

  col <- fig_cols()

  test_names <- names(me_dat$ma$obs)
  test_labels <- lapply(test_names, function(x) {
      tmp <- strsplit(toupper(x), split=".", fixed=TRUE);
      list(f=tmp[[1]][1], s=tmp[[1]][2])})
  names(test_labels) <- test_names

  par(mfrow=c(2,6))
  par(mar=c(4,1,1,1))
  for (x in names(me_dat$ma$obs)) {
    profiles.plot(me_dat$ma$sim[x], col.line=col[1],
                  opacity = 0.9, frame.plot=FALSE, yaxt="n",
                  xlab="", ylab="")
    abline(v=me_dat$ma$obs[,x], lty=2, lwd=2, col=col[2])
  }

  par(mar=c(4.5,1,1,1))
  for (x in names(fa_dat$ma$obs)) {
      ## set xlims properly
      if (fa_dat$ma$obs[x] > max(fa_dat$ma$sim[x]) | fa_dat$ma$obs[x] < min(fa_dat$ma$sim[x])){
          xlim <- c(as.numeric(fa_dat$ma$obs[x]-0.05), max(fa_dat$ma$sim[x]))
      } else {
          xlim <- range(fa_dat$ma$sim[x])
      }

      ## get labels for x axis
      label_f <- test_labels[x][[1]]$f
      label_s <- test_labels[x][[1]]$s
      xlab <- bquote(italic(.(label_f))[.(label_s)])

    profiles.plot(fa_dat$ma$sim[x], col.line=col[3],
                  opacity = 0.9, frame.plot=FALSE, yaxt="n",
                  xlab=xlab,
                  ylab="", cex.lab=1.5,
                  xlim=xlim)
    abline(v=fa_dat$ma$obs[,x], lty=2, lwd=2, col=col[2])
      
  }
}
