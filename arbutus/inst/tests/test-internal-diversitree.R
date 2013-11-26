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
lik.eb <- make.eb(phy, states)

fit.bm <- find.mle(lik.bm, .1)
fit.ou <- find.mle(lik.ou, c(.1, .1, mean(states)))
fit.eb <- find.mle(lik.eb, c(.1, 0))

set.seed(1)
samples.bm <- mcmc(lik.bm, coef(fit.bm), 100, w=0.1, print.every=0)
set.seed(1)
samples.ou <- mcmc(lik.ou, coef(fit.ou), 30,
                   w=c(0.1, 10, 10),
                   lower=c(0, 0, -10),
                   upper=c(Inf, 100, 10),
                   print.every=0)
set.seed(1)
samples.eb <- mcmc(lik.eb, coef(fit.eb), 30,
                   w=c(.2, 10), lower=c(0, -Inf), upper=c(Inf, 0),
                   print.every=0)

## Testing internal code; need some extra imports:
model.type <- arbutus:::model.type
model.data <- arbutus:::model.data
model.pars <- arbutus:::model.pars
model.info <- arbutus:::model.info
model.type.fit.mle <- arbutus:::model.type.fit.mle
model.data.fit.mle <- arbutus:::model.data.fit.mle
model.pars.fit.mle <- arbutus:::model.pars.fit.mle
model.info.fit.mle <- arbutus:::model.info.fit.mle
model.type.mcmcsamples <- arbutus:::model.type.mcmcsamples
model.data.mcmcsamples <- arbutus:::model.data.mcmcsamples
model.pars.mcmcsamples <- arbutus:::model.pars.mcmcsamples
model.info.mcmcsamples <- arbutus:::model.info.mcmcsamples

test_that("Model types are correct", {
  expect_that(model.type(fit.bm), is_identical_to("BM"))
  expect_that(model.type(fit.ou), is_identical_to("OU"))
  expect_that(model.type(fit.eb), is_identical_to("EB"))
  expect_that(model.type(samples.bm, lik.bm), is_identical_to("BM"))
  expect_that(model.type(samples.ou, lik.ou), is_identical_to("OU"))
  expect_that(model.type(samples.eb, lik.eb), is_identical_to("EB"))
})

test_that("Models return their source data", {
  cmp <- list(phy=phy, data=states[phy$tip.label])
  expect_that(model.data(fit.bm, lik.bm), equals(cmp))
  expect_that(model.data(fit.ou, lik.ou), equals(cmp))

  expect_that(model.data(samples.bm, lik.bm), equals(cmp))
  expect_that(model.data(samples.ou, lik.ou), equals(cmp))
})

test_that("Giving incorrect likelihood function fails", {
  states.false <- states * 1.1
  lik.bm.false <- make.bm(phy, states.false)
  lik.ou.false <- make.ou(phy, states.false)
  lik.eb.false <- make.eb(phy, states.false)

  expect_that(model.data(fit.bm, lik.bm.false), throws_error())
  expect_that(model.data(fit.ou, lik.ou.false), throws_error())
  expect_that(model.data(fit.eb, lik.eb.false), throws_error())

  ## Ommitting the likelihood function is an error
  expect_that(model.data(fit.bm), throws_error())
  expect_that(model.data(fit.ou), throws_error())
  expect_that(model.data(fit.eb), throws_error())

  ## But if we skip the checks, this should extract new false data:
  cmp.false <- list(phy=phy, data=states.false[phy$tip.label])
  expect_that(model.data(fit.bm, lik.bm.false, check=FALSE),
              equals(cmp.false))
  expect_that(model.data(fit.ou, lik.ou.false, check=FALSE),
              equals(cmp.false))
  expect_that(model.data(fit.eb, lik.eb.false, check=FALSE),
              equals(cmp.false))

  ## MCMC:
  expect_that(model.data(samples.bm, lik.bm.false), throws_error())
  expect_that(model.data(samples.ou, lik.ou.false), throws_error())
  expect_that(model.data(samples.eb, lik.eb.false), throws_error())

  ## Ommitting the likelihood function is an error
  expect_that(model.data(samples.bm), throws_error())
  expect_that(model.data(samples.ou), throws_error())
  expect_that(model.data(samples.eb), throws_error())

  ## But if we skip the checks, this should extract new false data:
  cmp.false <- list(phy=phy, data=states.false[phy$tip.label])
  expect_that(model.data(samples.bm, lik.bm.false, check=FALSE),
              equals(cmp.false))
  expect_that(model.data(samples.ou, lik.ou.false, check=FALSE),
              equals(cmp.false))
  expect_that(model.data(samples.eb, lik.eb.false, check=FALSE),
              equals(cmp.false))
})

test_that("Processed coefficient names are as expected", {
  expect_that(names(model.pars(fit.bm, lik.bm)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(model.pars(fit.ou, lik.ou)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(model.pars(fit.eb, lik.eb)),
              is_identical_to(arbutus:::parnames.eb()))

  expect_that(model.pars(fit.bm, lik.bm)$SE, equals(0))
  expect_that(model.pars(fit.ou, lik.ou)$SE, equals(0))
  expect_that(model.pars(fit.eb, lik.eb)$SE, equals(0))

  ## MCMC:
  expect_that(names(model.pars(samples.bm, lik.bm)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(model.pars(samples.ou, lik.ou)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(model.pars(samples.eb, lik.eb)),
              is_identical_to(arbutus:::parnames.eb()))

  expect_that(model.pars(samples.bm, lik.bm)$SE,
              equals(rep(0, nrow(samples.bm))))
  expect_that(model.pars(samples.ou, lik.ou)$SE,
              equals(rep(0, nrow(samples.ou))))
  expect_that(model.pars(samples.eb, lik.eb)$SE,
              equals(rep(0, nrow(samples.eb))))
})

test_that("Overall processed object looks legit", {
  obj.bm <- model.info(fit.bm, lik.bm)
  obj.ou <- model.info(fit.ou, lik.ou)
  obj.eb <- model.info(fit.eb, lik.eb)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))

  obj.bm <- model.info(samples.bm, lik.bm)
  obj.ou <- model.info(samples.ou, lik.ou)
  obj.eb <- model.info(samples.eb, lik.eb)

  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))
})

## Constrained functions
test_that("Constrained functions work", {
  lik.ou.c <- constrain(lik.ou, theta ~ 0)
  set.seed(1)
  fit.ou.c <- find.mle(lik.ou.c, c(.1, .1))
  samples.ou.c <- mcmc(lik.ou.c, coef(fit.ou.c), 30,
                       w=c(0.1, 10),
                       lower=c(0, 0),
                       upper=c(Inf, 100),
                       print.every=0)

  ## Can get the type out:
  expect_that(model.type(fit.ou.c), is_identical_to("OU"))
  expect_that(model.type(samples.ou.c, lik.ou.c),
              is_identical_to("OU"))

  ## Underlying data:
  cmp <- list(phy=phy, data=states[phy$tip.label])
  expect_that(model.data(fit.ou.c, lik.ou.c), equals(cmp))
  expect_that(model.data(samples.ou.c, lik.ou.c), equals(cmp))

  ## Parameters -- these need to include all parameters:
  p <- model.pars(fit.ou.c, lik.ou.c)
  expect_that(names(p), is_identical_to(arbutus:::parnames.ou()))
  expect_that(p$SE, equals(0))
  expect_that(p$z0, equals(0)) # z0 is theta, translated

  p <- model.pars(samples.ou.c, lik.ou.c)
  expect_that(names(p), is_identical_to(arbutus:::parnames.ou()))
  expect_that(p$SE, equals(rep(0, nrow(samples.ou))))
  expect_that(p$z0, equals(rep(0, nrow(samples.ou))))
})
