## Script for evaluating model adequacy of angiosperm clades
## Fitting models using maximum likelihood
source("R/modelfit-helper-fxns.R")

# This function does all the ML analysis.  It takes a data set (which
# contains a tree, character states, and standard errors for the
# character states), fits a model using ML and then assesses the
# adequacy of that model.
#
# Most of the ugliness is putting (hopefully) sensible bounds on the
# parameters to avoid ridges in likelihood space that can happen with
# OU and EB.
#
# Note that the entire adequacy part of the function is the line
#
#      ma <- phy.model.check(fit)  
model.ad.ml <- function(data, model) {
  model <- match.arg(model, c("BM", "OU", "EB"))
  ## Extract components from the pre-prepared data object.
  phy    <- data$phy
  states <- data$states
  SE     <- data$SE

  make.lik <- switch(model, BM=make.bm, OU=make.ou, EB=make.eb)
  lik <- make.lik(phy, states, SE,
                  control=list(method="pruning", backend="C"))

  # Start at REML estimate of sigsq for all models, using arbutus'
  # function to estiate this.
  s2 <- sigsq.est(as.unit.tree(phy, data=drop(states)))

  # ML bounds and starting points:
  lower <- 0
  upper <- Inf
  if (model == "BM") {
    start <- s2
  } else if (model == "OU") {
    upper <- bounds.ou(phy, states)
    start <- c(s2, 0.05)
  } else if (model == "EB") {
    lower <- c(0, -1)
    start <- c(s2, -0.1)
  }
  fit <- find.mle(lik, x.init=start)

  # Assess adequacy of all models
  ma <- phy.model.check(fit)

  info       <- data[c("taxa", "rank", "trait")]
  info$size  <- Ntip(phy)
  info$age   <- max(branching.times(phy))
  info$model <- model
  info$aic   <- AIC(fit)
  info$type  <- "ml"

  list(info=info, ma=ma)
}

dir.create(path.ml(), FALSE)

files <- dir(path.data())
ok <- lapply(files, run.model.ad, "ml", regenerate=TRUE)

# Next, process the output to create a file with summarised results.
fits <- lapply(dir(path.ml(), full.names=TRUE), readRDS)
results <- combine(fits)
write.csv(results, "output/results-ml.csv", row.names=FALSE)
