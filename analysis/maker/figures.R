fig_model_support_aic <- function(ml) {
  aic <- ml[aic_names]
  colnames(aic) <- c("BM", "OU", "EB")
  ## add dummy variable
  dd <- cbind(Subclade=as.character(seq_len(nrow(aic))), aic)
  ord <- dd[order(dd[,"OU"], dd[,"BM"], decreasing = TRUE), "Subclade"]
  ## reorient data frame for plotting
  df <- suppressMessages(melt(dd))
  colnames(df)[2:3] <- c("model", "weight")
  df$Subclade <- factor(df$Subclade, ord)

  ## create geom_bar plot
  .e <- environment()
  p <- ggplot(df, aes(factor(Subclade), weight, fill=model, order=model),
              environment=.e)
  p <- p + geom_bar(stat="identity", position="stack",width=1)
  p <- p + scale_y_continuous(name="AIC weight")
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
  p
}
