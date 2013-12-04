## For now, only BM models are implemented in diversitree.
model.type.fit.mle.pgls <- function(fit, ...) {
  "BM"
}

model.data.fit.mle.pgls <- function(fit, ..., check=TRUE) {
  require(diversitree)
  cache <- diversitree:::get.cache(get.likelihood(fit))
  list(phy=cache$info$phy, data=drop(resid(fit)))
}

model.pars.fit.mle.pgls <- function(fit, ..., check=TRUE) {
  ## TODO: Will only work with BM models; update once diversitree can
  ## do more than this.
  if (model.type(fit) != "BM")
    stop("Hard coded to deal with BM models only")
  s2 <- coef(fit, full=TRUE)[["s2"]]
  as.list(c(sigsq=s2, SE=0, z0=NA))
}
