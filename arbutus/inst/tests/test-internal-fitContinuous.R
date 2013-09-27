## Testing internal functions and testing functions in other packages
## are both considered to be poor form.  So this is something that we
## might re-address before release, but it will be good for now.
source("helper-arbutus.R")

## Testing internal code; need some extra imports:
modeltype <- arbutus:::modeltype
modeldata <- arbutus:::modeldata
modelpars <- arbutus:::modelpars
modelinfo <- arbutus:::modelinfo
modeltype.gfit <- arbutus:::modeltype.gfit
modeldata.gfit <- arbutus:::modeldata.gfit
modelpars.gfit <- arbutus:::modelpars.gfit
modelinfo.gfit <- arbutus:::modelinfo.gfit

context("fitContinuous (internals)")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[,"wingL"]
control <- list(niter=10)

set.seed(1)
fit.bm <- fitContinuous(phy, states, SE=NA, control=control, model="BM")
fit.ou <- suppressWarnings(fitContinuous(phy, states, SE=NA,
                                         control=control, model="OU"))
fit.eb <- fitContinuous(phy, states, SE=NA, control=control, model="EB")

## And without SE:
fit.bm.noSE <- fitContinuous(phy, states, control=control, model="BM")
fit.ou.noSE <- suppressWarnings(fitContinuous(phy, states,
                                              control=control, model="OU"))
fit.eb.noSE <- fitContinuous(phy, states, control=control, model="EB")

test_that("Model types are correct", {
  expect_that(modeltype(fit.bm), is_identical_to("BM"))
  expect_that(modeltype(fit.ou), is_identical_to("OU"))
  expect_that(modeltype(fit.eb), is_identical_to("EB"))

  expect_that(modeltype(fit.bm.noSE), is_identical_to("BM"))
  expect_that(modeltype(fit.ou.noSE), is_identical_to("OU"))
  expect_that(modeltype(fit.eb.noSE), is_identical_to("EB"))
})

test_that("Models return their source data", {
  cmp <- list(phy=phy, data=states)
  expect_that(modeldata(fit.bm), equals(cmp))
  expect_that(modeldata(fit.ou), equals(cmp))
  expect_that(modeldata(fit.eb), equals(cmp))

  expect_that(modeldata(fit.bm.noSE), equals(cmp))
  expect_that(modeldata(fit.ou.noSE), equals(cmp))
  expect_that(modeldata(fit.eb.noSE), equals(cmp))
})

test_that("Coefficient names are as expected", {
  pars.bm <- c("sigsq", "SE", "z0")
  pars.ou <- c("alpha", "sigsq", "SE", "z0")
  pars.eb <- c("a", "sigsq", "SE", "z0")
  
  expect_that(names(coef(fit.bm)), is_identical_to(pars.bm))
  expect_that(names(coef(fit.ou)), is_identical_to(pars.ou))
  expect_that(names(coef(fit.eb)), is_identical_to(pars.eb))

  expect_that(names(coef(fit.bm.noSE)),
              is_identical_to(setdiff(pars.bm, "SE")))
  expect_that(names(coef(fit.ou.noSE)),
              is_identical_to(setdiff(pars.ou, "SE")))
  expect_that(names(coef(fit.eb.noSE)),
              is_identical_to(setdiff(pars.eb, "SE")))
})

test_that("Processed coefficient names are as expected", {
  expect_that(names(modelpars(fit.bm)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(modelpars(fit.ou)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(modelpars(fit.eb)),
              is_identical_to(arbutus:::parnames.eb()))

  expect_that(modelpars(fit.bm)$SE, is_greater_than(0))
  expect_that(modelpars(fit.ou)$SE, is_greater_than(0))
  expect_that(modelpars(fit.eb)$SE, is_greater_than(0))
  
  expect_that(names(modelpars(fit.bm.noSE)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(modelpars(fit.ou.noSE)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(modelpars(fit.eb.noSE)),
              is_identical_to(arbutus:::parnames.eb()))

  expect_that(modelpars(fit.bm.noSE)$SE, equals(0))
  expect_that(modelpars(fit.ou.noSE)$SE, equals(0))
  expect_that(modelpars(fit.eb.noSE)$SE, equals(0))
})

test_that("Overall processed object looks legit", {
  obj.bm <- modelinfo(fit.bm)
  obj.ou <- modelinfo(fit.ou)
  obj.eb <- modelinfo(fit.eb)

  obj.bm.noSE <- modelinfo(fit.bm.noSE)
  obj.ou.noSE <- modelinfo(fit.ou.noSE)
  obj.eb.noSE <- modelinfo(fit.eb.noSE)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))

  expect_that(names(obj.bm.noSE), is_identical_to(obj.names))
  expect_that(names(obj.ou.noSE), is_identical_to(obj.names))
  expect_that(names(obj.eb.noSE), is_identical_to(obj.names))
})

test_that("Nonsense cases should fail", {
  expect_that(modelinfo(NULL), throws_error())
  expect_that(modelinfo(unclass(fit.bm)), throws_error())
})
