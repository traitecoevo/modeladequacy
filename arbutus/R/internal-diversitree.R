## See internal-fitContinuous.R for general information about what is
## going on here.

## TODO: check issue #9.

model.type.fit.mle <- function(fit, ...) {
  toupper(sub("^fit\\.mle\\.", "", class(fit)[[1]]))
}

model.data.fit.mle <- function(fit, ...) {
  require(diversitree)
  cache <- diversitree:::get.cache(get.likelihood(fit))
  list(phy=cache$info$phy, data=cache$states)
}

model.pars.fit.mle <- function(fit, ...) {
  require(diversitree)
  model <- model.type(fit)
  pars <- as.list(coef(fit, full=TRUE))
  pars <- diversitree.to.arbutus.recast(pars, model)
  pars <- diversitree.to.arbutus.se(pars, get.likelihood(fit))
  pars
}

#' @method model.info fit.mle
#' @S3method model.info fit.mle
model.info.fit.mle <- function(fit, ...) {
  m <- list(data=model.data(fit),
            pars=model.pars(fit),
            type=model.type(fit))
  class(m) <- "fitC"
  m
}

## Similar, but for mcmcsamples:
model.type.mcmcsamples <- function(fit, ...) {
  require(diversitree)
  lik <- get.likelihood(fit)
  if (is.constrained(lik))
    toupper(class(get.likelihood(lik))[[1]])
  else
    toupper(class(lik)[[1]])
}

model.data.mcmcsamples <- function(fit, ...) {
  require(diversitree)
  cache <- diversitree:::get.cache(get.likelihood(fit))
  list(phy=cache$info$phy, data=cache$states)
}

model.pars.mcmcsamples <- function(fit, burnin=NA, thin=NA, sample=NA, ...) {
  require(diversitree)
  model <- model.type(fit)
  lik <- get.likelihood(fit)
  pars <- as.data.frame(coef(fit, burnin=burnin, thin=thin,
                             sample=sample, full=TRUE, lik=lik))
  pars <- diversitree.to.arbutus.recast(pars, model)
  pars <- diversitree.to.arbutus.se(pars, lik)
  pars
}

#' @method model.info mcmcsamples
#' @S3method model.info mcmcsamples
model.info.mcmcsamples <- function(fit, burnin=NA, thin=NA, sample=NA, ...) {
  m <- list(data=model.data(fit),
            pars=model.pars(fit, burnin, thin, sample),
            type=model.type(fit))
  class(m) <- c("multiFitC", "fitC")
  m
}

diversitree.to.arbutus.recast <- function(pars, model) {
  if (!(model %in% c("BM", "OU", "EB")))
    stop(sprintf("Cannot extract from models of model %s", model))

  tr <- c(sigsq="s2", z0="theta")
  idx <- match(names(pars), tr)
  i <- !is.na(idx)
  names(pars)[i] <- names(tr)[idx][i]

  if (model == "OU") {
    pars <- pars[c("alpha", "sigsq", "z0")]
  } else if (model == "EB") {
    pars <- pars[c("a", "sigsq")]
    pars$z0 <- NA
  } else if (model == "BM") {
    ## At this point we don't have access to the root state.  Given
    ## that different root treatments are possible in diversitree,
    ## this is not actually very straightforward, either.
    pars$z0 <- NA
  }

  pars
}

diversitree.to.arbutus.se <- function(pars, lik) {
  require(diversitree)
  ## Pull the standard error out.  Also ugly.
  cache <- diversitree:::get.cache(lik)
  sd <- cache$states.sd
  if (length(unique(sd)) == 1)
    pars$SE <- sd[[1]]
  else
    stop("Variable length states.sd -- cannot deal with this (yet)")
  pars
}
