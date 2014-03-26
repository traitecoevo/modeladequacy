## Functions used in analyses
## But not part of main analyses

## Load in all packages used 
library(arbutus)
library(diversitree)
library(multicore)

# TODO: Find and replace manual uses of angio-data with this.
path.data <- function() {
  "output/angio-data"
}
path.ml <- function() {
  "output/results-ml"
}
path.bayes <- function() {
  "output/results-bayes"
}

# Little wrapper function that loads data, fits a model (ML or MCMC),
# checks the model adequacy, and then saves the results in a file.  If
# the output file exists, this is skipped.
run.model.ad <- function(filename, type=c("ml", "bayes"),
                         regenerate=FALSE, verbose=TRUE) {
  type <- match.arg(type)
  if (type == "ml") {
    path.out <- path.ml()
    model.ad <- model.ad.ml
  } else {
    path.out <- path.bayes()
    model.ad <- model.ad.bayes
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
    res <- lapply(models, function(m) model.ad(dat, m))
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
    # P values for the various statistics:
    pv <- pval.summ.stats(x$ma)
    names(pv) <- sprintf("%s.%s.%s", type, tolower(x$info$model), names(pv))

    # Mean distance
    #
    # TODO: This fails sometimes (computationally singular) but I
    # think I'm out of date with arbutus.  For now, just setting this
    # to use try(), and setting to NaN on error so we can track it
    # down later.
    mv <- try(mean(mv.summ.stats(x$ma)), silent=TRUE)
    if (inherits(mv, "try-error")) {
      mv <- NaN
    }
    names(mv) <- sprintf("mv.%s.%s", type, tolower(x$info$model))
    
    c(list(dic=x$info[[col.ic]]),
      as.list(pv),
      as.list(mv))
  }
  
  info <- res$BM$info[c("taxa", "rank", "trait", "size", "age")]

  weights <- as.list(ic.weights(sapply(res, function(x) x$info[[col.ic]])))
  names(weights) <- sprintf("%sw.%s", col.ic, tolower(names(weights)))

  tmp <- lapply(res, process1)
  stats <- unlist(tmp, recursive=FALSE, use.names=FALSE)
  names(stats) <- unlist(lapply(tmp, names))
  
  as.data.frame(c(info, weights, stats), stringsAsFactors=FALSE)
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
