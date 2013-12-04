library(diversitree)
library(nlme)
library(testthat)

context("GLS (diversitree) (internals)")

model.type <- arbutus:::model.type
model.data <- arbutus:::model.data
model.pars <- arbutus:::model.pars
model.info <- arbutus:::model.info

model.type.fit.mle <- arbutus:::model.type.fit.mle
model.data.fit.mle <- arbutus:::model.data.fit.mle
model.pars.fit.mle <- arbutus:::model.pars.fit.mle
model.info.fit.mle <- arbutus:::model.info.fit.mle

model.type.fit.mle.pgls <- arbutus:::model.type.fit.mle.pgls
model.data.fit.mle.pgls <- arbutus:::model.data.fit.mle.pgls
model.pars.fit.mle.pgls <- arbutus:::model.pars.fit.mle.pgls
# model.info.fit.mle.pgls <- arbutus:::model.info.fit.mle.pgls

set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
tx <- sim.character(phy, 1)
ty <- sim.character(phy, 1) + 3
data <- data.frame(x=tx, y=ty, row.names=names(tx))

fit.gls.bm.ml   <- gls(y ~ x, data, corBrownian(1, phy), method="ML")
fit.gls.bm.reml <- gls(y ~ x, data, corBrownian(1, phy), method="REML")

lik.pgls.bm.vcv <- make.pgls(phy, y ~ x, data,
                             control=list(method="vcv"))
lik.pgls.bm.con <- make.pgls(phy, y ~ x, data,
                             control=list(method="contrasts"))

s2.ml <- arbutus:::estimate.sigma2.gls(fit.gls.bm.ml)
p.ml <- c(coef(fit.gls.bm.ml), s2=s2.ml)


fit.pgls.bm.vcv <- find.mle(lik.pgls.bm.vcv, c(0, 0, 1))
fit.pgls.bm.con <- find.mle(lik.pgls.bm.con, c(0, 0, 1))

test_that("The fits are sane", {
  expect_that(fit.pgls.bm.vcv, equals(fit.pgls.bm.con))
})

test_that("Model types are correct", {
  expect_that(model.type(fit.pgls.bm.vcv), is_identical_to("BM"))
  expect_that(model.type(fit.pgls.bm.con), is_identical_to("BM"))
})

test_that("Models return their source data", {
  drop.attributes <- function(x) {
    attributes(x) <- attributes(x)["names"]
    x
  }
  cmp.bm <- list(phy=phy, data=drop.attributes(resid(fit.gls.bm.ml)))
  expect_that(model.data(fit.pgls.bm.vcv), equals(cmp.bm, tolerance=1e-5))
  expect_that(model.data(fit.pgls.bm.con), equals(cmp.bm, tolerance=1e-5))

  ## Update the fits to use the same values as gls() and we'll get
  ## better agreement.
  fit.pgls.bm.vcv$par <- p.ml
  fit.pgls.bm.con$par <- p.ml
  expect_that(model.data(fit.pgls.bm.vcv), equals(cmp.bm, tolerance=1e-13))
  expect_that(model.data(fit.pgls.bm.con), equals(cmp.bm, tolerance=1e-13))
})

test_that("Model parameters are correct", {
  cmp <- as.list(c(sigsq=coef(fit.pgls.bm.vcv)[["s2"]], SE=0, z0=NA))
  expect_that(model.pars(fit.pgls.bm.vcv), equals(cmp))
  expect_that(model.pars(fit.pgls.bm.con), equals(cmp))
})

test_that("Overall processed object looks legit", {
  obj.bm.vcv <- model.info(fit.pgls.bm.vcv)
  obj.bm.con <- model.info(fit.pgls.bm.con)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm.vcv), is_identical_to(obj.names))
  expect_that(names(obj.bm.con), is_identical_to(obj.names))
})
