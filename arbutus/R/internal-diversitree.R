## See internal-fitContinuous.R for general information about what is
## going on here.

## So, despite what we spoke about earlier, diversitree does not
## include the likelihood function within the fitted object -- see
## issue #9.  This means that many of these require both the fitted
## object and the likelihood function, but will by default check that
## the likelihood function could plausibly been used to fit the model.
model.type.fit.mle <- function(fit, ...) {
  toupper(sub("^fit\\.mle\\.", "", class(fit)[[1]]))
}

model.data.fit.mle <- function(fit, lik, ..., check=TRUE) {
   ## require diversitree
   require(diversitree)
   if (check)
    check.diversitree.lik.fit(fit, lik, error=TRUE, ...)
  cache <- diversitree:::get.cache(lik)
  list(phy=cache$info$phy, data=cache$states)
}

model.pars.fit.mle <- function(fit, lik, ..., check=TRUE) {
  if (check)
    check.diversitree.lik.fit(fit, lik, error=TRUE, ...)
  model <- model.type(fit)
  pars <- as.list(coef(fit, full=TRUE))
  pars <- diversitree.to.arbutus.recast(pars, model)
  pars <- diversitree.to.arbutus.se(pars, lik)
  pars
}

#' @method model.info fit.mle
#' @S3method model.info fit.mle
model.info.fit.mle <- function(fit, lik, ..., check=TRUE) {
  m <- list(data=model.data(fit, lik, check=check),
            pars=model.pars(fit, lik, check=FALSE), # no need to check twice
            type=model.type(fit))
  class(m) <- "fitC"
  m
}

check.diversitree.lik.fit <- function(fit, lik, error=TRUE, ...) {
  ll <- lik(coef(fit), ...)
  ok <- isTRUE(all.equal(ll, fit$lnLik, check.attributes=FALSE))
  if (error && !ok)
    stop("Likelihood function and fitted object do not agree")
  invisible(ok)
}

## Similar, but for mcmcsamples:
model.type.mcmcsamples <- function(fit, lik, ..., check=TRUE) {
  if (check)
    check.diversitree.lik.fit.mcmcsamples(fit, lik, error=TRUE, ...)
  if (is.constrained(lik))
    toupper(class(attr(lik, "func"))[[1]])
  else
    toupper(class(lik)[[1]])
}

model.data.mcmcsamples <- function(fit, lik, ..., check=TRUE) {
   ## require diversitree
   require(diversitree)
   if (check)
    check.diversitree.lik.fit.mcmcsamples(fit, lik, error=TRUE, ...)
  cache <- diversitree:::get.cache(lik)
  list(phy=cache$info$phy, data=cache$states)
}

model.pars.mcmcsamples <- function(fit, lik, ..., check=TRUE) {
  if (check)
    check.diversitree.lik.fit.mcmcsamples(fit, lik, error=TRUE, ...)
  model <- model.type(fit, lik, check=FALSE) # don't check again.

  pars <- as.data.frame(coef(fit, full=TRUE, lik=lik))
  pars <- diversitree.to.arbutus.recast(pars, model)
  pars <- diversitree.to.arbutus.se(pars, lik)
  pars
}

#' @method model.info mcmcsamples
#' @S3method model.info mcmcsamples
model.info.mcmcsamples <- function(fit, lik, ..., check=TRUE) {
  m <- list(data=model.data(fit, lik, check=check),
            pars=model.pars(fit, lik, check=FALSE), # no need to check twice
            type=model.type(fit, lik, check=FALSE))
  class(m) <- c("multiFitC", "fitC")
  m
}

## Bit of a hack.
check.diversitree.lik.fit.mcmcsamples <- function(fit, lik, error=TRUE,
                                                  ..., prior=NULL) {
  idx <- 1
  p <- coef(fit)[idx,]
  ll <- lik(p, ...)
  if (!is.null(prior))
    ll <- ll + prior(p)
  ok <- isTRUE(all.equal(ll, fit$p[idx], check.attributes=FALSE))
  if (error && !ok)
    stop("Likelihood function and fitted object do not agree")
  invisible(ok)
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
  ## require diversitree
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
