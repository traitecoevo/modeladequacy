## Functions used in analyses
## But not part of main analyses

## Load in all packages used 
library(arbutus)
library(diversitree)
library(parallel)

# This function does all the ML or Bayesian analysis, depending on the
# 'type' argument.  It takes a data set (which contains a tree,
# character states, and standard errors for the character states),
# fits a model using ML and then assesses the adequacy of that model.
#
# Most of the ugliness is putting (hopefully) sensible bounds on the
# parameters to avoid ridges in likelihood space that can happen with
# OU and EB.
#
# Note that the entire adequacy part of the function is the line
#      ma <- phy.model.check(fit)
# or
#     ma <- phy.model.check(samples)
model.ad <- function(data, model, type, seed=1) {
  model <- match.arg(model, c("BM", "OU", "EB"))
  type  <- match.arg(type,  c("ml", "bayes"))
  ## Extract components from the pre-prepared data object.
  phy    <- data$phy
  states <- drop(data$states)
  SE     <- data$SE

  # Make the analyses recomputable by using the same seed each time:
  set.seed(seed)

  make.lik <- switch(model, BM=make.bm, OU=make.ou, EB=make.eb)
  lik <- make.lik(phy, states, SE,
                  control=list(method="pruning", backend="C"))

  # Start at REML estimate of sigsq for all models, using arbutus'
  # function to estiate this.
  s2 <- sigsq.est(as.unit.tree(phy, data=drop(states)))

  # ML bounds, starting points and priors (priors only used in the
  # Bayesian analysis)
  lower <- 0
  upper <- Inf
  if (model == "BM") {
    start <- s2
    prior <- make.prior.bm(s2.lower=0, s2.upper=2)
  } else if (model == "OU") {
    upper <- bounds.ou(phy, states)
    start <- c(s2, 0.05)
    prior <- make.prior.ou(s2.lower=0, s2.upper=2,
                           ln.mean=log(0.5), ln.sd=log(1.5))
  } else if (model == "EB") {
    lower <- c(0, -1)
    start <- c(s2, -0.1)
    prior <- make.prior.eb(s2.lower=0, s2.upper=2,
                           a.lower=-1, a.upper=0)
  }

  if (type == "ml") {
    fit <- find.mle(lik, x.init=start)
    ic  <- AIC(fit)
    ic.name <- "aic"

    # Assess adequacy of all models
    ma <- phy.model.check(fit)
  } else if (type == "bayes") {
    # Some general parameters
    pilot  <- 100
    burnin <- 1000
    nsteps <- 10000

    # Run short chain to obtain appropriate step size for MCMC
    tmp <- mcmc(lik, x.init=start, nsteps=pilot, prior=prior, w=1,
                print.every=0)

    w <- diff(apply(coef(tmp), 2, range))

    # Full chain:
    samples <- mcmc(lik, x.init=start, nsteps=nsteps, prior=prior, w=w,
                    print.every=0)

    ic <- dic.mcmcsamples(samples, burnin=burnin)
    ic.name <- "dic"

    # Assess adequacy of all models    
    ma <- phy.model.check(samples, burnin=burnin, sample=1000)
  }

  # For all cases, the general information about three we want is the same:
  info       <- data[c("taxa", "rank", "trait")]
  info$size  <- Ntip(phy)
  info$age   <- max(branching.times(phy))
  info$model <- model
  info$type  <- type
  info[[ic.name]] <- ic

  # Returning information about the tree, the assessment of model
  # adequacy and the fit (either the ML object or the MCMC samples).
  list(info=info, ma=ma, fit=if(type == "ml") fit else samples)
}

# Little wrapper function that loads data, fits a model (ML or MCMC),
# checks the model adequacy, and then saves the results in a file.  If
# the output file exists, this is skipped.
run.model.ad <- function(filename, type=c("ml", "bayes"),
                         regenerate=FALSE, verbose=TRUE) {
  type <- match.arg(type)
  if (type == "ml") {
    path.out <- path.ml()
  } else {
    path.out <- path.bayes()
  }

  filename.out <- file.path(path.out, filename)
  if (file.exists(filename.out) && !regenerate) {
    invisible(FALSE)
  } else {
    models <- c("BM", "OU", "EB")
    dat <- readRDS(file.path(path.data(), filename))
    if (verbose) {
      message(sprintf("%s: %s", type, sub(".rds$", "", filename)))
    }
    res <- lapply(models, function(m) model.ad(dat, m, type))
    names(res) <- models
    saveRDS(res, filename.out)
    invisible(TRUE)
  }
}


combine <- function(fits) {
  do.call(rbind, lapply(fits, process))
}

process <- function(res) {
  type <- res$BM$info$type
  col.ic <- if (type == "ml") "aic" else "dic"

  process1 <- function(x) {
    model <- tolower(x$info$model)
    # P values for the various statistics:
    pv <- pval.summ.stats(x$ma)
    names(pv) <- sprintf("%s.%s.%s", names(pv), type, model)

    # Mean distance
    #
    # NOTE: This fails sometimes (computationally singular) but I
    # think I'm out of date with arbutus.  For now, just setting this
    # to use try(), and setting to NaN on error so we can track it
    # down later.
    mv <- try(mean(mv.summ.stats(x$ma)), silent=TRUE)
    if (inherits(mv, "try-error")) {
      mv <- NaN
    }
    names(mv) <- sprintf("mv.%s.%s", type, model)

    # Mean diagnostic
    #
    # This is not strictly necessary but was used to test whether
    # m.pic was consistently under or overestimated This is the only
    # summary statistic which we expected to be systematically biased
    # in one direction
    md <- sum(x$ma$summ.stats.sim$m.pic < x$ma$summ.stats.obs$m.pic)
    names(md) <- sprintf("mean.diag.%s", model)

    as.list(c(pv, mv, md))
  }
  
  info <- res$BM$info[c("taxa", "rank", "trait", "size", "age")]

  ic <- sapply(res, function(x) x$info[[col.ic]])
  icw <- ic.weights(ic)
  names(ic) <- sprintf("%s.%s", col.ic, tolower(names(ic)))
  names(icw) <- sprintf("%sw.%s", col.ic, tolower(names(icw)))

  stats <- unlist.with.names(lapply(res, process1))
  
  as.data.frame(c(info, ic, icw, stats), stringsAsFactors=FALSE)
}

## Compute weights for AIC or DIC
## takes a named vector of AIC/DIC values
ic.weights <- function(x){
  d.x <- x - min(x)
  tmp <- exp(-0.5 * d.x)
  tmp / sum(tmp)
}

## Compute deviance information criterion from mcmcsamples
dic.mcmcsamples <- function(x, burnin=0){
  if (!inherits(x, "mcmcsamples"))
    stop("this function is only designed for diversitree's mcmcsamples object")

  lik <- attr(x, "func")
  p <- coef(x, burnin=burnin)
  dev <- -2 * apply(p, 1, lik)
  
  ## estimate effective number of parameters
  dbar <- mean(dev)

  ## deviance of posterior means:
  post.means <- colMeans(p)
  ## evaluate deviance at the mean posterior estimate
  dhat <- -2 * lik(post.means)

  pd <- dbar - dhat

  ## calculate dic
  dic <- dbar + pd
  unname(dic)
}


## function for making prior for bm
## just a wrapper for make.prior.uniform
make.prior.bm <- function(s2.lower, s2.upper) {
  make.prior.uniform(lower=s2.lower, upper=s2.upper)
}



## define function to create a lognormal prior
## to be used in make.prior.ou
make.prior.lognormal <- function(ln.mean, ln.sd) {
  function(x) sum(dlnorm(x, meanlog=ln.mean, sdlog=ln.sd, log=TRUE))
}



## function for making ou prior
## s2 is given a uniform prior
## alpha is given a lognormal prior
make.prior.ou <- function(s2.lower, s2.upper, ln.mean, ln.sd){
  p.s2 <- make.prior.uniform(lower=s2.lower, upper=s2.upper)
  p.al <- make.prior.lognormal(ln.mean=ln.mean, ln.sd=ln.sd)

  function(pars) {
    p.s2(pars[[1]]) + p.al(pars[[2]])
  }
}



## function for making eb prior
## s2 is given a uniform prior
## a is givne a uniform prior
## the prior on s2 should be different than that of a
## s2 > 0 and a < 0
make.prior.eb <- function(s2.lower, s2.upper, a.lower, a.upper){
  p.s2 <- make.prior.uniform(lower=s2.lower, upper=s2.upper)
  p.a  <- make.prior.uniform(lower=a.lower, upper=a.upper)

  function(pars) {
    p.s2(pars[[1]]) + p.a(pars[[2]])
  }
}


## create function for setting OU bounds for ML analyses
## this is the only one where bounds to be necessary
bounds.ou <- function(phy, states){
  ## find the shortest terminal branch
  sht <- min(phy$edge.length[phy$edge[,2] <= Ntip(phy)])
  ## upper bound for alpha: greater than shortest branch
  ## upper bound for for sigsq: set to 2
  c(2, log(2)/sht)
}

unlist.with.names <- function(x) {
  ret <- unlist(x, recursive=FALSE, use.names=FALSE)
  names(ret) <- unlist(lapply(x, names))
  ret
}
