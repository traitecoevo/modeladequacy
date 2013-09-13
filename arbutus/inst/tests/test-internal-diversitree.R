## Testing internal functions and testing functions in other packages
## are both considered to be poor form.  So this is something that we
## might re-address before release, but it will be good for now.
source("helper-arbutus.R")

context("diversitree (internals)")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[,"wingL"]

lik.bm <- make.bm(phy, states)
lik.ou <- make.ou(phy, states)

fit.bm <- find.mle(lik.bm, .1)
fit.ou <- find.mle(lik.ou, c(.1, .1, mean(states)))

## Testing internal code; need some extra imports:
modeltype <- arbutus:::modeltype
modeldata <- arbutus:::modeldata
modelpars <- arbutus:::modelpars
modelinfo <- arbutus:::modelinfo
modeltype.fit.mle <- arbutus:::modeltype.fit.mle
modeldata.fit.mle <- arbutus:::modeldata.fit.mle
modelpars.fit.mle <- arbutus:::modelpars.fit.mle
modelinfo.fit.mle <- arbutus:::modelinfo.fit.mle

test_that("Model types are correct", {
  modeltype(fit.bm, is_identical_to("BM"))
  modeltype(fit.ou, is_identical_to("OU"))
})

test_that("Models return their source data", {
  cmp <- list(phy=phy, data=states[phy$tip.label])
  expect_that(modeldata(fit.bm, lik.bm), equals(cmp))
  expect_that(modeldata(fit.ou, lik.ou), equals(cmp))
})

test_that("Giving incorrect likelihood function fails", {
  states.false <- states * 1.1
  lik.bm.false <- make.bm(phy, states.false)
  lik.ou.false <- make.ou(phy, states.false)

  expect_that(modeldata(fit.bm, lik.bm.false), throws_error())
  expect_that(modeldata(fit.ou, lik.ou.false), throws_error())

  ## Ommitting the likelihood function is an error
  expect_that(modeldata(fit.bm), throws_error())
  expect_that(modeldata(fit.ou), throws_error())

  ## But if we skip the checks, this should extract new false data:
  cmp.false <- list(phy=phy, data=states.false[phy$tip.label])
  expect_that(modeldata(fit.bm, lik.bm.false, check=FALSE),
              equals(cmp.false))
  expect_that(modeldata(fit.ou, lik.ou.false, check=FALSE),
              equals(cmp.false))
})

test_that("Processed coefficient names are as expected", {
  expect_that(names(modelpars(fit.bm, lik.bm)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(modelpars(fit.ou, lik.ou)),
              is_identical_to(arbutus:::parnames.ou()))

  expect_that(modelpars(fit.bm, lik.bm)$SE, equals(0))
  expect_that(modelpars(fit.ou, lik.ou)$SE, equals(0))
})

test_that("Overall processed object looks legit", {
  obj.bm <- modelinfo(fit.bm, lik.bm)
  obj.ou <- modelinfo(fit.ou, lik.ou)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
})
