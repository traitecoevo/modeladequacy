source("helper-arbutus.R")

context("GLS (internals)")

model.type <- arbutus:::model.type
model.data <- arbutus:::model.data
model.pars <- arbutus:::model.pars
model.info <- arbutus:::model.info

model.type.gls <- arbutus:::model.type.gls
model.data.gls <- arbutus:::model.data.gls
model.pars.gls <- arbutus:::model.pars.gls
model.info.gls <- arbutus:::model.info.gls

set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
tx <- sim.character(phy, 1)
ty <- sim.character(phy, 1) + 3
data <- data.frame(x=tx, y=ty, row.names=names(tx))

## Note that the EB model won't fit without the fixed=TRUE bit added.
fit.gls.bm <- gls(y ~ x, data, corBrownian(1, phy))
fit.gls.ou <- gls(y ~ x, data, corMartins(1, phy))
fit.gls.eb <- gls(y ~ x, data, corBlomberg(1, phy, fixed=TRUE))
fit.gls.la <- gls(y ~ x, data, corPagel(1, phy))

fit.gls.bm.ml <- gls(y ~ x, data, corBrownian(1, phy), method="ML")
fit.gls.ou.ml <- gls(y ~ x, data, corMartins(1, phy), method="ML")
fit.gls.eb.ml <- gls(y ~ x, data, corBlomberg(1, phy, fixed=TRUE), method="ML")
fit.gls.la.ml <- gls(y ~ x, data, corPagel(1, phy), method="ML")

test_that("Model types are correct", {
  expect_that(model.type(fit.gls.bm), is_identical_to("BM"))
  expect_that(model.type(fit.gls.ou), is_identical_to("OU"))
  expect_that(model.type(fit.gls.eb), is_identical_to("EB"))
  expect_that(model.type(fit.gls.la), is_identical_to("lambda"))
})

test_that("Models return their source data", {
  drop.attributes <- function(x) {
    attributes(x) <- attributes(x)["names"]
    x
  }
  cmp.bm <- list(phy=phy, data=drop.attributes(resid(fit.gls.bm)))
  cmp.ou <- list(phy=phy, data=drop.attributes(resid(fit.gls.ou)))
  cmp.eb <- list(phy=phy, data=drop.attributes(resid(fit.gls.eb)))
  cmp.la <- list(phy=phy, data=drop.attributes(resid(fit.gls.la)))

  expect_that(model.data(fit.gls.bm), equals(cmp.bm))
  expect_that(model.data(fit.gls.ou), equals(cmp.ou))
  expect_that(model.data(fit.gls.eb), equals(cmp.eb))
  expect_that(model.data(fit.gls.la), equals(cmp.la))
})

## Checking that the s2 values are correct is actually quite
## difficult.  A decent check here will have to wait.  In the interim,
## here is Matt's original function, so at least it's a separate
## implementation to check against.
get.sigma2.gls <- function(fit) {
  model <- model.type(fit)
  phy <- model.data(fit)$phy
  rr <- as.numeric(resid(fit))
  names(rr) <- names(resid(fit))
   
  if (model == "BM"){
    sigsq <- sigsq.est(as.unit.tree(phy, data=rr))
    pars <- list(sigsq=sigsq, z0=NA, SE=0)
  } else if (model == "OU"){
    alpha <- as.numeric(fit$modelStruct$corStruct)
    tmp <- list(sigsq=1, alpha=alpha, SE=0)
    rescalephy <- model.phylo.ou(phy, pars=tmp)
    sigsq <- sigsq.est(as.unit.tree(rescalephy, data=rr)) 
  } else if (model == "EB"){
    a <- as.numeric(fit$modelStruct$corStruct)
    tmp <- list(sigsq=1, a=a, SE=0)
    rescalephy <- model.phylo.eb(phy, pars=tmp)
    sigsq <- sigsq.est(as.unit.tree(rescalephy, data=rr)) 
  } else if (model == "lambda") {
    lambda <- as.numeric(fit$modelStruct$corStruct)
    tmp <- list(sigsq=1, lambda=lambda, SE=0)
    rescalephy <- model.phylo.lambda(phy, pars=tmp)
    sigsq <- sigsq.est(as.unit.tree(rescalephy, data=rr))  
  }

  if (fit$method == "ML")
    sigsq <- sigsq * (fit$dims$N - fit$dims$p) / fit$dims$N

  sigsq
}
environment(get.sigma2.gls) <- environment(arbutus:::sigsq.est)
estimate.sigma2.gls <- arbutus:::estimate.sigma2.gls

test_that("Diffusion parameter is correct (REML)", {
  expect_that(estimate.sigma2.gls(fit.gls.bm),
              equals(get.sigma2.gls(fit.gls.bm)))
  expect_that(estimate.sigma2.gls(fit.gls.ou),
              equals(get.sigma2.gls(fit.gls.ou)))
  expect_that(estimate.sigma2.gls(fit.gls.eb),
              equals(get.sigma2.gls(fit.gls.eb)))
  expect_that(estimate.sigma2.gls(fit.gls.la),
              equals(get.sigma2.gls(fit.gls.la)))
})

test_that("Diffusion parameter is correct (ML)", {
  expect_that(estimate.sigma2.gls(fit.gls.bm.ml),
              equals(get.sigma2.gls(fit.gls.bm.ml)))
  expect_that(estimate.sigma2.gls(fit.gls.ou.ml),
              equals(get.sigma2.gls(fit.gls.ou.ml)))
  expect_that(estimate.sigma2.gls(fit.gls.eb.ml),
              equals(get.sigma2.gls(fit.gls.eb.ml)))
  expect_that(estimate.sigma2.gls(fit.gls.la.ml),
              equals(get.sigma2.gls(fit.gls.la.ml)))

  ## This is the only relationship where this is guaranteed, I
  ## believe.
  expect_that(estimate.sigma2.gls(fit.gls.bm.ml),
              is_less_than(estimate.sigma2.gls(fit.gls.bm)))
})

test_that("Coefficient names are as expected", {
  pars.bm <- c("sigsq", "SE", "z0")
  pars.ou <- c("alpha", "sigsq", "SE", "z0")
  pars.eb <- c("a", "sigsq", "SE", "z0")
  pars.la <- c("lambda", "sigsq", "SE", "z0")

  expect_that(names(model.pars(fit.gls.bm)), is_identical_to(pars.bm))
  expect_that(names(model.pars(fit.gls.ou)), is_identical_to(pars.ou))
  expect_that(names(model.pars(fit.gls.eb)), is_identical_to(pars.eb))
  expect_that(names(model.pars(fit.gls.la)), is_identical_to(pars.la))
})

test_that("Overall processed object looks legit", {
  obj.bm <- model.info(fit.gls.bm)
  obj.ou <- model.info(fit.gls.ou)
  obj.eb <- model.info(fit.gls.eb)
  obj.la <- model.info(fit.gls.la)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))
  expect_that(names(obj.la), is_identical_to(obj.names))

  expect_that(class(obj.bm), is_identical_to("fitC"))
  expect_that(class(obj.ou), is_identical_to("fitC"))
  expect_that(class(obj.eb), is_identical_to("fitC"))
  expect_that(class(obj.la), is_identical_to("fitC"))
})
