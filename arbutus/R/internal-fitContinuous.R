## Internal functions for pulling apart fitContinuous objects.  These
## need to be kept in close alignment with geiger, because they use
## undocumented (and therefore potentially unstable) features of
## geiger.

## Some of these could be merged into geiger, but we're going to need
## similar features for every model type and people won't do this for
## us.  So instead, we'll write some unit tests that check that things
## are behaving as expected.  Then we should at least get fairly early
## warning when things break.  Later on, we can look at rolling some
## interface into the different packages if people will let us.

## Note that we depend on version 2.0 of geiger, and previous versions
## *will not work*, here.

## Note even though these aren't actually S3 generics (yet) the
## function names are written as though they are.
modeltype.gfit <- function(fit, ...) {
  ret <- attr(fit$lik, "model")
  if (is.null(ret) || !is.character(ret) || length(ret) != 1)
    stop("Failed to extract model type from geiger fit")
  ret
}

modeldata.gfit <- function(fit, ...) {
  e <- environment(fit$lik)
  list(phy=get("phy", e), data=get("dat", e))
}

modelpars.gfit <- function(fit, ...) {
  model <- modeltype(fit)
    pars <- as.list(coef(fit))
    if (!("SE" %in% names(pars)))
      pars$SE <- 0
    pars[c(setdiff(names(pars), "SE"), "SE")]
}

modelinfo.gfit <- function(fit, ...) {
  m <- list(data=modeldata(fit),
            pars=modelpars(fit),
            type=modeltype(fit))
  class(m) <- "model.fitC"
  m
}
