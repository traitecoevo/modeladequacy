source("helper-arbutus.R")
library(caper)

context("Unit tree")

pic.var <- function(...)
  pic(..., var.contrasts=TRUE)

## First, work with a unit tree built directly from a tree and state
## data; this is generally not what is done.  In fact, the name here
## is a bit surprising, because this does not really make a unit tree;
## that is done by the rescaling functions.
##
## TODO: Think about whether `as.unit.tree.phylo` is correctly named.
data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- dat$phy
states <- dat$dat[phy$tip.label,"wingL"]

test_that("Unit tree construction works", {
  phy.unit <- as.unit.tree(phy, states)

  expect_that(phy.unit, is_a("unit.tree"))
  expect_that(names(phy.unit),
              is_identical_to(c("phy", "data", "pics")))
  expect_that(phy.unit$phy, is_a("phylo"))
  expect_that(phy.unit$data, is_a("numeric"))
  expect_that(phy.unit$pics, is_a("matrix"))

  ## Compare with exactly what a unit tree has:
  expect_that(phy.unit$phy, is_identical_to(phy))
  expect_that(phy.unit$data, is_identical_to(states))
  expect_that(phy.unit$pics,
              is_identical_to(pic.var(states, phy)))

  ## Check that the is.unit.tree function works
  expect_that(is.unit.tree(phy.unit), is_true())
})

## Simple tests; we require sanitised tree and data coming in;
## everything within as.unit.tree.phylo can rely on that.
test_that("Invalid inputs fail appropriately", {
  expect_that(as.unit.tree(phy, dat$dat[,"wingL"]), throws_error())
  expect_that(as.unit.tree(NULL, states), throws_error())
  expect_that(as.unit.tree(phy, NULL), throws_error())
  expect_that(as.unit.tree(phy, cbind(states)), throws_error())
})

## Test all the different models that might be used with a unit tree.
## This is quite a lot of different cases, and making sure that they
## are all tested is a trick.  One of the things we're really looking
## for here is that argument passing works as expected.  Of course,
## most methods don't really take arguments, so that's OK.
##
## Most other options should have been done by lower level tests.
test_that("Unit tree construction from BM/fitContinuous works", {
  fit.bm <- fitContinuousQuiet(phy=phy, dat=states, model="BM",
                               control=list(niter=5))
  phy.unit <- as.unit.tree(fit.bm)
  expect_that(phy.unit, is_a("unit.tree"))

  ## Manual rescaling using Geiger function
  cmp <- rescale(phy, "BM", coef(fit.bm)[[1]])
  expect_that(phy.unit$phy, equals(cmp))
  expect_that(phy.unit$data, equals(states))
  expect_that(phy.unit$pics, equals(pic.var(states, cmp)))
})

test_that("Unit tree construction from BM/fitContinuous with SE works", {
  se <- 0.01
  fit.bm <- fitContinuousQuiet(phy=phy, dat=states, model="BM", SE=se,
                               control=list(niter=5))
  phy.unit <- as.unit.tree(fit.bm)
  expect_that(phy.unit, is_a("unit.tree"))

  ## Manual rescaling using Geiger function
  cmp <- rescale.se(phy, "BM", coef(fit.bm)[[1]], SE=se)
  expect_that(phy.unit$phy, equals(cmp))
  expect_that(phy.unit$data, equals(states))
  expect_that(phy.unit$pics, equals(pic.var(states, cmp)))
})

test_that("Unit tree construction from BM/diversitree/ML works", {
  lik.bm <- make.bm(phy, states)
  fit.bm <- find.mle(lik.bm, .1)
  phy.unit <- as.unit.tree(fit.bm)
  expect_that(phy.unit, is_a("unit.tree"))
  
  ## Manual rescaling using Geiger function
  cmp <- rescale(phy, "BM", coef(fit.bm)[[1]])
  expect_that(phy.unit$phy, equals(cmp))
  expect_that(phy.unit$data, equals(states))
  expect_that(phy.unit$pics, equals(pic.var(states, cmp)))
})

test_that("Unit tree construction from BM/diversitree/ML with SE works", {
  se <- 0.01
  lik.bm <- make.bm(phy, states, states.sd=se)
  fit.bm <- find.mle(lik.bm, .1)
  phy.unit <- as.unit.tree(fit.bm)
  expect_that(phy.unit, is_a("unit.tree"))
  
  ## Manual rescaling using Geiger function
  cmp <- rescale.se(phy, "BM", coef(fit.bm)[[1]], SE=se)
  expect_that(phy.unit$phy, equals(cmp))
  expect_that(phy.unit$data, equals(states))
  expect_that(phy.unit$pics, equals(pic.var(states, cmp)))
})

test_that("Unit tree construction from BM/diversitree/MCMC works", {
  lik.bm <- make.bm(phy, states)
  fit.bm <- find.mle(lik.bm, .1)
  set.seed(1)
  samples.bm <- mcmc(lik.bm, coef(fit.bm), 100, w=0.1, print.every=0)

  phy.unit <- as.unit.tree(samples.bm)
  expect_that(phy.unit, is_a("multiPhylo"))
  expect_that(phy.unit[[1]], is_a("unit.tree"))
  expect_that(length(phy.unit), equals(nrow(samples.bm)))

  ## Quick check:
  idx <- 5
  fit.bm$par <- coef(samples.bm)[idx,]
  fit.bm$lnLik <- samples.bm$p[idx]
  cmp <- as.unit.tree(fit.bm)
  expect_that(phy.unit[[idx]], is_identical_to(cmp))

  ## Check burn-in and samples work
  nb <- 10
  sample <- 20
  expect_that(length(as.unit.tree(samples.bm, burnin=nb)),
              equals(nrow(coef(samples.bm, burnin=nb))))
  expect_that(length(as.unit.tree(samples.bm, sample=sample)),
              equals(nrow(coef(samples.bm, sample=sample))))
})

test_that("Unit tree construction from BM/diversitree/MCMC with SE works", {
  se <- 0.01
  lik.bm <- make.bm(phy, states, states.sd=se)
  fit.bm <- find.mle(lik.bm, .1)
  set.seed(1)
  samples.bm <- mcmc(lik.bm, coef(fit.bm), 100, w=0.1, print.every=0)

  phy.unit <- as.unit.tree(samples.bm)
  expect_that(phy.unit, is_a("multiPhylo"))
  expect_that(phy.unit[[1]], is_a("unit.tree"))
  expect_that(length(phy.unit), equals(nrow(samples.bm)))

  ## Quick check:
  idx <- 5
  fit.bm$par <- coef(samples.bm)[idx,]
  fit.bm$lnLik <- samples.bm$p[idx]
  cmp <- as.unit.tree(fit.bm)
  expect_that(phy.unit[[idx]], is_identical_to(cmp))

  ## Check burn-in and samples work
  nb <- 10
  sample <- 20
  expect_that(length(as.unit.tree(samples.bm, burnin=nb)),
              equals(nrow(coef(samples.bm, burnin=nb))))
  expect_that(length(as.unit.tree(samples.bm, sample=sample)),
              equals(nrow(coef(samples.bm, sample=sample))))
})

## Now, we quickly blast through some of the other tests.
##
## A full list would be:
##   fitContinuous: (BM, OU, EB, lambda, kappa, delta, white) x (with/without SE)
##   diversitree: (BM, OU, EB) x (with/without SE) x (ML, MCMC)
##   caper (pgls): lambda, kappa, gamma
##   gls: BM, OU, EB, lambda
##   phylolm: BM, OU, EB, lambda, kappa, delta
##   diversitree pgls: BM x (ML, MCMC)
##
## which is a lot (7 * 2 + 3 * 2 * 2 + 3 + 4 + 6 + 2 = 41 models)
##
## I consider it to be an open question what the best way of dealing
## with this is, given that all of these models should be fairly well
## tested by their tests to `model.info`.
##
## The tests below are really checking the tree rescaling more than
## anything, and that is separately tested in test-rescale.R
test_that("OU tree rescaling worked (fitContinuous)", {
    fit.ou <- fitContinuousQuiet(phy=phy, dat=states, model="OU",
                                 control=list(niter=10))
    phy.unit <- as.unit.tree(fit.ou)

    ## Manual rescaling using Geiger function
    cmp <- rescale(phy, "OU", alpha=coef(fit.ou)[[1]], sigsq=coef(fit.ou)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("OU tree rescaling worked (diversitree)", {
    lik.ou <- make.ou(phy, states)
    fit.ou <- find.mle(lik.ou, x.init=c(0.1,0.1))

    phy.unit <- as.unit.tree(fit.ou)

    ## Manual rescaling using Geiger function
    cmp <- rescale(phy, "OU", alpha=coef(fit.ou)[[2]], sigsq=coef(fit.ou)[[1]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("OU tree rescaling worked (diversitree, mcmc)", {
  lik.ou <- make.ou(phy, states)
  fit.ou <- find.mle(lik.ou, x.init=c(0.1,0.1))
  set.seed(1)
  samples.ou <- mcmc(lik.ou, coef(fit.ou), 30,
                     w=c(0.1, 10),
                     lower=c(0, 0),
                     upper=c(Inf, 100),
                     print.every=0)

  phy.unit <- as.unit.tree(samples.ou)
  expect_that(phy.unit, is_a("multiPhylo"))
  expect_that(phy.unit[[1]], is_a("unit.tree"))

  ## Quick check:
  idx <- 5
  fit.ou$par <- coef(samples.ou)[idx,]
  fit.ou$lnLik <- samples.ou$p[idx]
  cmp <- as.unit.tree(fit.ou)
  expect_that(phy.unit[[idx]], is_identical_to(cmp))
})

test_that("EB tree rescaling worked (fitContinuous)", {
    fit.eb <- fitContinuousQuiet(phy=phy, dat=states, model="EB",
                                 control=list(niter=10))
    phy.unit <- as.unit.tree(fit.eb)

    ## Manual rescaling using Geiger functions
    cmp <- rescale(phy, "EB", a=coef(fit.eb)[[1]], sigsq=coef(fit.eb)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("lambda tree rescaling worked (fitContinuous)", {
    fit.lamb <- fitContinuousQuiet(phy=phy, dat=states, model="lambda",
                                   control=list(niter=10))
    phy.unit <- as.unit.tree(fit.lamb)

    ## Manual rescaling using Geiger functions
    cmp <- rescale(phy, "lambda", lambda=coef(fit.lamb)[[1]], sigsq=coef(fit.lamb)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("kappa tree rescaling worked (fitContinuous)", {
    fit.k <- fitContinuousQuiet(phy=phy, dat = states, model="kappa",
                                control = list(niter=10))

    phy.unit <- as.unit.tree(fit.k)

    ## Manual rescaling using Geiger functions
    cmp <- rescale(phy, "kappa", kappa=coef(fit.k)[[1]], sigsq=coef(fit.k)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("delta tree rescaling worked (fitContinuous)", {
    fit.d <- fitContinuousQuiet(phy=phy, dat=states, model="delta",
                                control = list(niter=10))

    phy.unit <- as.unit.tree(fit.d)

    ## Manual rescaling using Geiger function
    ## Using internal fxn as rescale.phy does a post-hoc adjustment
    ## We are not using adjustment here
    foo <- geiger:::.delta.phylo(phy)
    cmp <- foo(delta=coef(fit.d)[[1]], sigsq=coef(fit.d)[[2]], rescale=FALSE)

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("white noise rescaling worked (fitContinuous)", {
    fit.w <- fitContinuousQuiet(phy=phy, dat=states, model="white",
                                control = list(niter=10))

    phy.unit <- as.unit.tree(fit.w)

    ## Manual rescaling using Geiger function
    cmp <- rescale(phy, "white", sigsq=coef(fit.w)[[1]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

## The PGLS code is a bit different, especially the diversitree mcmc
## code, so we're going to test that separately and more thoroughly
## than the code above.

## PGLS
set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
tx <- sim.character(phy, 1)
ty <- sim.character(phy, 1) + 3
data <- data.frame(x=tx, y=ty, row.names=names(tx))

## nlme::gls:
fit.gls.bm.ml   <- gls(y ~ x, data, corBrownian(1, phy), method="ML")
fit.gls.bm.reml <- gls(y ~ x, data, corBrownian(1, phy), method="REML")

## caper::pgls:
cdata <- comparative.data(phy, cbind(data, sp=rownames(data)), 'sp')
fit.caper <- pgls(y ~ x, cdata)

## Diversitree, ML
lik.pgls.bm.vcv <- make.pgls(phy, y ~ x, data,
                             control=list(method="vcv"))
lik.pgls.bm.con <- make.pgls(phy, y ~ x, data,
                             control=list(method="contrasts"))

s2.ml   <- arbutus:::estimate.sigma2.gls(fit.gls.bm.ml)
s2.reml <- arbutus:::estimate.sigma2.gls(fit.gls.bm.reml)
p.ml <- c(coef(fit.gls.bm.ml), s2=s2.ml)

fit.pgls.bm.vcv <- find.mle(lik.pgls.bm.vcv, c(0, 0, 1))
fit.pgls.bm.con <- find.mle(lik.pgls.bm.con, c(0, 0, 1))

## Draw some mcmc samples, straight from the likelihood.
set.seed(1)
samples.pgls <- mcmc(lik.pgls.bm.con, p.ml, 100, w=1, print.every=0)

drop.attributes <- function(x) {
  at <- attributes(x)
  attributes(x) <- at["names"]
  x
}

## Now, start testing.
test_that("Unit tree construction from gls works (ML)", {
  phy.unit <- as.unit.tree(fit.gls.bm.ml)
  expect_that(phy.unit, is_a("unit.tree"))

  cmp <- rescale(phy, "BM", s2.ml)
  dat <- drop.attributes(resid(fit.gls.bm.ml))
  
  expect_that(phy.unit$phy,  equals(cmp))
  expect_that(phy.unit$data, equals(dat))
  expect_that(phy.unit$pics, equals(pic.var(dat, cmp)))
})

test_that("Unit tree construction from gls works (REML)", {
  phy.unit <- as.unit.tree(fit.gls.bm.reml)
  expect_that(phy.unit, is_a("unit.tree"))

  cmp <- rescale(phy, "BM", s2.reml)
  dat <- drop.attributes(resid(fit.gls.bm.reml))
  
  expect_that(phy.unit$phy,  equals(cmp))
  expect_that(phy.unit$data, equals(dat))
  expect_that(phy.unit$pics, equals(pic.var(dat, cmp)))
})

test_that("Unit tree construction from caper pgls works", {
  phy.unit <- as.unit.tree(fit.caper)
  expect_that(phy.unit, is_a("unit.tree"))

  cmp <- rescale(phy, "BM", s2.ml)
  dat <- drop(resid(fit.caper))
  
  expect_that(phy.unit$phy,  equals(cmp))
  expect_that(phy.unit$data, equals(dat))
  expect_that(phy.unit$pics, equals(pic.var(dat, cmp)))
})

## TODO: The ML s2 estimate is incorrect - missing the -1 term off k.
test_that("Unit tree construction from diversitre pgls (vcv) works", {
  phy.unit <- as.unit.tree(fit.pgls.bm.vcv)
  expect_that(phy.unit, is_a("unit.tree"))

  s2 <- coef(fit.pgls.bm.vcv)[["s2"]]
  cmp <- rescale(phy, "BM", s2)
  dat <- drop(resid(fit.pgls.bm.vcv))

  expect_that(phy.unit$phy,  equals(cmp))
  expect_that(phy.unit$data, equals(dat))
  expect_that(phy.unit$pics, equals(pic.var(dat, cmp)))
})


## TODO: The ML s2 estimate is incorrect - missing the -1 term off k.
test_that("Unit tree construction from diversitre pgls (contrasts) works", {
  phy.unit <- as.unit.tree(fit.pgls.bm.con)
  expect_that(phy.unit, is_a("unit.tree"))

  s2 <- coef(fit.pgls.bm.con)[["s2"]]
  cmp <- rescale(phy, "BM", s2)
  dat <- drop(resid(fit.pgls.bm.con))

  expect_that(phy.unit$phy,  equals(cmp))
  expect_that(phy.unit$data, equals(dat))
  expect_that(phy.unit$pics, equals(pic.var(dat, cmp)))
})

test_that("Unit tree construction from diversitre pgls (mcmc) works", {
  phy.unit <- as.unit.tree(samples.pgls)
  expect_that(phy.unit, is_a("multiPhylo"))
  expect_that(phy.unit[[1]], is_a("unit.tree"))
  expect_that(length(phy.unit), equals(nrow(samples.pgls)))

  ## Quick check:
  idx <- 5
  s2 <- samples.pgls$s2[[idx]]
  cmp <- rescale(phy, "BM", s2)
  dat <- resid(samples.pgls)[,idx]

  expect_that(phy.unit[[idx]]$phy,  equals(cmp))
  expect_that(phy.unit[[idx]]$data, equals(dat))
  expect_that(phy.unit[[idx]]$pics, equals(pic.var(dat, cmp)))

  nb <- 10
  sample <- 20
  expect_that(length(as.unit.tree(samples.pgls, burnin=nb)),
              equals(nrow(coef(samples.pgls, burnin=nb))))
  expect_that(length(as.unit.tree(samples.pgls, sample=sample)),
              equals(nrow(coef(samples.pgls, sample=sample))))
})
