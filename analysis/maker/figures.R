fig_cols <- function() {
  c("#F46D43", "#3288BD", "#CDCD00")
}

fig_model_support_ic <- function(fits) {
  type <- if ("aic.bm" %in% names(fits)) "aic" else "dic"
  col <- fig_cols()
  ylab <- paste(toupper(type), "weight")

  col_names <- paste0(type, "w.", c("bm", "ou", "eb"))
  ic <- fits[col_names]

  colnames(ic) <- c("BM", "OU", "EB")
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
  df$trait <- factor(df$trait,
                     levels=c("sla", "seed_mass", "leaf_n"),
                     labels=c("SLA", "SeedMass", "LeafN"))

  df$variable <- factor(df$variable,
                        levels=pvalue_names_arbutus(),
                        labels=c("italic(M[SIG])", "italic(C[VAR])", "italic(S[VAR])", "italic(S[ASR])", "italic(S[HGT])", "italic(D[CDF])"))
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
