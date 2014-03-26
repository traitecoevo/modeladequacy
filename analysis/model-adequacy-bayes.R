## Script for evaluating model adequacy of angiosperm clades
## Fitting models using maximum likelihood

## Read in all the helper files
source("R/modelfit-helper-fxns.R")

model.ad.bayes <- function(data, model) {
  # NOTE: A lot of overlap here between this and the model.ad.ml()
  # function, but we're leaving it in for didactic purposes.
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

  if (model == "BM") {
    prior <- make.prior.bm(s2.lower=0, s2.upper=2)
    start <- s2
  } else if (model == "OU") {
    prior <- make.prior.ou(s2.lower=0, s2.upper=2,
                           ln.mean=log(0.5), ln.sd=log(1.5))
    start <- c(s2, 0.05)
  } else if (model == "EB") {
    prior <- make.prior.eb(s2.lower=0, s2.upper=2,
                           a.lower=-1, a.upper=0)
    start <- c(s2, -0.1)
  }

  # Some general parameters
  pilot  <- 100
  burnin <- 1000
  nsteps <- 10000

  # Make the analyses recomputable by using the same seed each time:
  set.seed(1)
  
  # Run short chain to obtain appropriate step size for MCMC
  tmp <- mcmc(lik, x.init=start, nsteps=pilot, prior=prior, w=1,
              print.every=0)

  w <- diff(apply(coef(tmp), 2, range))

  # Full chain:
  samples <- mcmc(lik, x.init=start, nsteps=nsteps, prior=prior, w=w,
                  print.every=0)

  dic <- dic.mcmcsamples(samples, burnin=burnin)
  ma <- phy.model.check(samples, burnin=burnin, sample=1000)

  info       <- data[c("taxa", "rank", "trait")]
  info$size  <- Ntip(phy)
  info$age   <- max(branching.times(phy))
  info$model <- model
  info$dic   <- dic.mcmcsamples(samples, burnin=burnin)
  info$type  <- "bayes"

  list(info=info, ma=ma)
}

dir.create(path.bayes(), FALSE)

files <- dir(path.data())
ok <- lapply(files, run.model.ad, "bayes", regenerate=TRUE)

# Next, process the output to create a file with summarised results.
fits <- lapply(dir(path.bayes(), full.names=TRUE), readRDS)
results <- combine(fits)
write.csv(results, "output/results-bayes.csv", row.names=FALSE)
