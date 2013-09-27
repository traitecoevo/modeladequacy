## Testing internal functions and testing functions in other packages
## are both considered to be poor form.  So this is something that we
## might re-address before release, but it will be good for now.
source("helper-arbutus.R")

## Testing internal code; need some extra imports:
model.type <- arbutus:::model.type
model.data <- arbutus:::model.data
model.pars <- arbutus:::model.pars
model.info <- arbutus:::model.info
model.type.gfit <- arbutus:::model.type.gfit
model.data.gfit <- arbutus:::model.data.gfit
model.pars.gfit <- arbutus:::model.pars.gfit
model.info.gfit <- arbutus:::model.info.gfit

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
fit.eb <- suppressWarnings(fitContinuous(phy, states, SE=NA,
                                         control=control, model="EB"))
fit.lam <- suppressWarnings(fitContinuous(phy, states, SE=NA,
                                          control=control, model="lambda"))
fit.del <- suppressWarnings(fitContinuous(phy, states, SE=NA,
                                          control=control, model="delta"))
fit.kap <- suppressWarnings(fitContinuous(phy, states, SE=NA,
                                          control=control, model="kappa"))
fit.wh <- suppressWarnings(fitContinuous(phy, states, SE=NA,
                                         control=control, model="white"))

## And without SE:
fit.bm.noSE <- fitContinuous(phy, states, control=control, model="BM")
fit.ou.noSE <- suppressWarnings(fitContinuous(phy, states,
                                              control=control, model="OU"))
fit.eb.noSE <- suppressWarnings(fitContinuous(phy, states,
                                         control=control, model="EB"))
fit.lam.noSE <- suppressWarnings(fitContinuous(phy, states,
                                          control=control, model="lambda"))
fit.del.noSE <- suppressWarnings(fitContinuous(phy, states,
                                          control=control, model="delta"))
fit.kap.noSE <- suppressWarnings(fitContinuous(phy, states,
                                          control=control, model="kappa"))
fit.wh.noSE <- suppressWarnings(fitContinuous(phy, states,
                                         control=control, model="white"))


test_that("Model types are correct", {
  expect_that(model.type(fit.bm), is_identical_to("BM"))
  expect_that(model.type(fit.ou), is_identical_to("OU"))
  expect_that(model.type(fit.eb), is_identical_to("EB"))
  expect_that(model.type(fit.lam), is_identical_to("lambda"))
  expect_that(model.type(fit.del), is_identical_to("delta"))
  expect_that(model.type(fit.kap), is_identical_to("kappa"))
  expect_that(model.type(fit.wh), is_identical_to("white"))

  expect_that(model.type(fit.bm.noSE), is_identical_to("BM"))
  expect_that(model.type(fit.ou.noSE), is_identical_to("OU"))
  expect_that(model.type(fit.eb.noSE), is_identical_to("EB"))
  expect_that(model.type(fit.lam.noSE), is_identical_to("lambda"))
  expect_that(model.type(fit.del.noSE), is_identical_to("delta"))
  expect_that(model.type(fit.kap.noSE), is_identical_to("kappa"))
  expect_that(model.type(fit.wh.noSE), is_identical_to("white"))
})

test_that("Models return their source data", {
  cmp <- list(phy=phy, data=states)
  expect_that(model.data(fit.bm), equals(cmp))
  expect_that(model.data(fit.ou), equals(cmp))
  expect_that(model.data(fit.eb), equals(cmp))
  expect_that(model.data(fit.lam), equals(cmp))
  expect_that(model.data(fit.del), equals(cmp))
  expect_that(model.data(fit.kap), equals(cmp))
  expect_that(model.data(fit.wh), equals(cmp))

  expect_that(model.data(fit.bm.noSE), equals(cmp))
  expect_that(model.data(fit.ou.noSE), equals(cmp))
  expect_that(model.data(fit.eb.noSE), equals(cmp))
  expect_that(model.data(fit.lam.noSE), equals(cmp))
  expect_that(model.data(fit.del.noSE), equals(cmp))
  expect_that(model.data(fit.kap.noSE), equals(cmp))
  expect_that(model.data(fit.wh.noSE), equals(cmp))

  
})

test_that("Coefficient names are as expected", {
  pars.bm <- c("sigsq", "SE", "z0")
  pars.ou <- c("alpha", "sigsq", "SE", "z0")
  pars.eb <- c("a", "sigsq", "SE", "z0")
  pars.lam <- c("lambda", "sigsq", "SE", "z0")
  pars.kap <- c("kappa", "sigsq", "SE", "z0")
  pars.del <- c("delta", "sigsq", "SE", "z0")
  pars.wh <- c("sigsq", "SE", "z0")
  
  expect_that(names(coef(fit.bm)), is_identical_to(pars.bm))
  expect_that(names(coef(fit.ou)), is_identical_to(pars.ou))
  expect_that(names(coef(fit.eb)), is_identical_to(pars.eb))
  expect_that(names(coef(fit.lam)), is_identical_to(pars.lam))
  expect_that(names(coef(fit.del)), is_identical_to(pars.del))
  expect_that(names(coef(fit.kap)), is_identical_to(pars.kap))
  expect_that(names(coef(fit.wh)), is_identical_to(pars.wh))

  expect_that(names(coef(fit.bm.noSE)),
              is_identical_to(setdiff(pars.bm, "SE")))
  expect_that(names(coef(fit.ou.noSE)),
              is_identical_to(setdiff(pars.ou, "SE")))
  expect_that(names(coef(fit.eb.noSE)),
              is_identical_to(setdiff(pars.eb, "SE")))
  expect_that(names(coef(fit.lam.noSE)),
              is_identical_to(setdiff(pars.lam, "SE")))
  expect_that(names(coef(fit.del.noSE)),
              is_identical_to(setdiff(pars.del, "SE")))
  expect_that(names(coef(fit.kap.noSE)),
              is_identical_to(setdiff(pars.kap, "SE")))
  expect_that(names(coef(fit.wh.noSE)),
              is_identical_to(setdiff(pars.wh, "SE")))  
})

test_that("Processed coefficient names are as expected", {
  expect_that(names(model.pars(fit.bm)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(model.pars(fit.ou)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(model.pars(fit.eb)),
              is_identical_to(arbutus:::parnames.eb()))
  expect_that(names(model.pars(fit.lam)),
              is_identical_to(arbutus:::parnames.lambda()))
  expect_that(names(model.pars(fit.del)),
              is_identical_to(arbutus:::parnames.delta()))
  expect_that(names(model.pars(fit.kap)),
              is_identical_to(arbutus:::parnames.kappa()))
  expect_that(names(model.pars(fit.wh)),
              is_identical_to(arbutus:::parnames.white()))
 

  expect_that(model.pars(fit.bm)$SE, is_greater_than(0))
  expect_that(model.pars(fit.ou)$SE, is_greater_than(0))
  expect_that(model.pars(fit.eb)$SE, is_greater_than(0))
  expect_that(model.pars(fit.lam)$SE, is_greater_than(0))
  expect_that(model.pars(fit.del)$SE, is_greater_than(0))
  expect_that(model.pars(fit.kap)$SE, is_greater_than(0))
  expect_that(model.pars(fit.wh)$SE, is_greater_than(0))
  
  expect_that(names(model.pars(fit.bm.noSE)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(model.pars(fit.ou.noSE)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(model.pars(fit.eb.noSE)),
              is_identical_to(arbutus:::parnames.eb()))
  expect_that(names(model.pars(fit.lam.noSE)),
              is_identical_to(arbutus:::parnames.lambda()))
  expect_that(names(model.pars(fit.del.noSE)),
              is_identical_to(arbutus:::parnames.delta()))
  expect_that(names(model.pars(fit.kap.noSE)),
              is_identical_to(arbutus:::parnames.kappa()))
  expect_that(names(model.pars(fit.wh.noSE)),
              is_identical_to(arbutus:::parnames.white()))
  

  expect_that(model.pars(fit.bm.noSE)$SE, equals(0))
  expect_that(model.pars(fit.ou.noSE)$SE, equals(0))
  expect_that(model.pars(fit.eb.noSE)$SE, equals(0))
  expect_that(model.pars(fit.lam.noSE)$SE, equals(0))
  expect_that(model.pars(fit.del.noSE)$SE, equals(0))
  expect_that(model.pars(fit.kap.noSE)$SE, equals(0))
  expect_that(model.pars(fit.wh.noSE)$SE, equals(0))
  
})

test_that("Overall processed object looks legit", {
  obj.bm <- model.info(fit.bm)
  obj.ou <- model.info(fit.ou)
  obj.eb <- model.info(fit.eb)
  obj.lam <- model.info(fit.lam)
  obj.del <- model.info(fit.del)
  obj.kap <- model.info(fit.kap)
  obj.wh <- model.info(fit.wh)

  obj.bm.noSE <- model.info(fit.bm.noSE)
  obj.ou.noSE <- model.info(fit.ou.noSE)
  obj.eb.noSE <- model.info(fit.eb.noSE)
  obj.lam.noSE <- model.info(fit.lam.noSE)
  obj.kap.noSE <- model.info(fit.kap.noSE)
  obj.del.noSE <- model.info(fit.del.noSE)
  obj.wh.noSE <- model.info(fit.wh.noSE)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))
  expect_that(names(obj.lam), is_identical_to(obj.names))
  expect_that(names(obj.del), is_identical_to(obj.names))
  expect_that(names(obj.kap), is_identical_to(obj.names))
  expect_that(names(obj.wh), is_identical_to(obj.names))

  expect_that(names(obj.bm.noSE), is_identical_to(obj.names))
  expect_that(names(obj.ou.noSE), is_identical_to(obj.names))
  expect_that(names(obj.eb.noSE), is_identical_to(obj.names))
  expect_that(names(obj.lam.noSE), is_identical_to(obj.names))
  expect_that(names(obj.del.noSE), is_identical_to(obj.names))
  expect_that(names(obj.kap.noSE), is_identical_to(obj.names))
  expect_that(names(obj.wh.noSE), is_identical_to(obj.names))
})

test_that("Nonsense cases should fail", {
  expect_that(modelinfo(NULL), throws_error())
  expect_that(modelinfo(unclass(fit.bm)), throws_error())
})
