## See internal-fitContinuous.R for general information about what is
## going on here.

## So, despite what we spoke about earlier, diversitree does not
## include the likelihood function within the fitted object -- see
## issue #9.  This means that many of these require both the fitted
## object and the likelihood function, but will by default check that
## the likelihood function could plausibly been used to fit the model.
modeltype.fit.mle <- function(fit, ...) {
  toupper(sub("^fit\\.mle\\.", "", class(fit)[[1]]))
}

modeldata.fit.mle <- function(fit, lik, ..., check=TRUE) {
  if (check)
    check.diversitree.lik.fit(fit, lik, error=TRUE)
  cache <- diversitree:::get.cache(lik)
  list(phy=cache$info$phy, data=cache$states)
}

modelpars.fit.mle <- function(fit, lik, ..., check=TRUE) {
  if (check)
    check.diversitree.lik.fit(fit, lik, error=TRUE)
  model <- modeltype(fit)
  if (model %in% c("BM", "OU")) {
    pars <- as.list(coef(fit))

    ## Translate parameter names (ugly)
    tr <- c(sigsq="s2", z0="theta")
    idx <- match(names(pars), tr)
    i <- !is.na(idx)
    names(pars)[i] <- names(tr)[idx][i]

    ## Reorder parameter names, and fill in missing ones.
    if (model == "OU")
      pars <- pars[c("alpha", "sigsq", "z0")]
    else if (model == "BM")
      ## At this point we don't have access to the root state.  Given
      ## that different root treatments are possible in diversitree,
      ## this is not actually very straightforward, either.
      pars$z0 <- NA

    ## Pull the standard error out.  Also ugly.
    cache <- diversitree:::get.cache(lik)
    sd <- cache$states.sd
    if (length(unique(sd)) == 1)
      pars$SE <- sd[[1]]
    else
      stop("Variable length states.sd -- cannot deal with this (yet)")

    pars
  } else {
    stop(sprintf("Cannot extract from models of model %s", model))
  }
}


## As this will use the rescaling fxns from fitContinuous, I have
## called this an object of "model.fitC" class
modelinfo.fit.mle <- function(fit, lik, ..., check=TRUE) {
  m <- list(data=modeldata(fit, lik, check=check),
            pars=modelpars(fit, lik, check=FALSE), # no need to check twice
            type=modeltype(fit))
  class(m) <- "model.fitC"
  m
}

check.diversitree.lik.fit <- function(fit, lik, error=TRUE, ...) {
  ll <- lik(coef(fit), ...)
  ok <- isTRUE(all.equal(ll, fit$lnLik, check.attributes=FALSE))
  if (error && !ok)
    stop("Likelihood function and fitted object do not agree")
  invisible(ok)
}
