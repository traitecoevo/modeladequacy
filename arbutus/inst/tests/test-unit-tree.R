source("helper-arbutus.R")

context("Unit tree")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[,"wingL"]

test_that("Unit tree construction works", {
  phy.unit <- as.unit.tree(phy, states)
  expect_that(names(phy.unit),
              is_identical_to(c("phy", "data", "pics")))
  expect_that(phy.unit$phy, is_a("phylo"))
  expect_that(phy.unit$data, is_a("matrix"))
  expect_that(phy.unit$pics, is_a("matrix"))

  expect_that(phy.unit$phy, is_identical_to(phy))
  expect_that(phy.unit$data,
              equals(cbind(states, deparse.level=0)))
  expect_that(phy.unit$pics,
              equals(pic(cbind(states), phy, var.contrasts=TRUE)))
})

test_that("BM tree rescaling worked (fitContinuous)", {
  fit.bm <- fitContinuous(phy=phy, dat=states, model="BM",
                          control=list(niter=10))
  phy.unit <- as.unit.tree(fit.bm)

  ## Manual rescaling using Geiger function
  cmp <- rescale(phy, "BM", coef(fit.bm)[[1]])

  expect_that(phy.unit$phy, equals(cmp))
  expect_that(phy.unit$phy$edge.length / phy$edge.length,
              equals(rep(coef(fit.bm)[[1]], length(phy$edge.length))))
  expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("BM tree rescaling worked (diversitree)", {
  lik.bm <- make.bm(phy, states)
  fit.bm <- find.mle(lik.bm, .1)
  phy.unit <- as.unit.tree(fit.bm, lik.bm)

  cmp <- rescale(phy, "BM", coef(fit.bm)[[1]])

  expect_that(phy.unit$phy, equals(cmp))
  expect_that(phy.unit$phy$edge.length / phy$edge.length,
              equals(rep(coef(fit.bm)[[1]], length(phy$edge.length))))
  expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("OU tree rescaling worked (fitContinuous)", {
    fit.ou <- suppressWarnings(fitContinuous(phy=phy, dat=states, model="OU",
                            control=list(niter=10)))
    phy.unit <- as.unit.tree(fit.ou)

    ## Manual rescaling using Geiger function
    cmp <- rescale(phy, "OU", alpha=coef(fit.ou)[[1]], sigsq=coef(fit.ou)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("OU tree rescaling worked (diversitree)", {
    lik.ou <- make.ou(phy, states)
    fit.ou <- find.mle(lik.ou, x.init=c(0.1,0.1, mean(states)))

    phy.unit <- as.unit.tree(fit.ou, lik.ou)

    ## Manual rescaling using Geiger function
    cmp <- rescale(phy, "OU", alpha=coef(fit.ou)[[2]], sigsq=coef(fit.ou)[[1]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("EB tree rescaling worked (fitContinuous)", {
    fit.eb <- suppressWarnings(fitContinuous(phy=phy, dat=states, model="EB",
                                             control=list(niter=10)))
    phy.unit <- as.unit.tree(fit.eb)

    ## Manual rescaling using Geiger functions
    cmp <- rescale(phy, "EB", a=coef(fit.eb)[[1]], sigsq=coef(fit.eb)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("lambda tree rescaling worked (fitContinuous)", {
    fit.lamb <- suppressWarnings(fitContinuous(phy=phy, dat=states, model="lambda",
                                               control=list(niter=10)))
    phy.unit <- as.unit.tree(fit.lamb)

    ## Manual rescaling using Geiger functions
    cmp <- rescale(phy, "lambda", lambda=coef(fit.lamb)[[1]], sigsq=coef(fit.lamb)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("kappa tree rescaling worked (fitContinuous)", {
    fit.k <- suppressWarnings(fitContinuous(phy=phy, dat = states, model="kappa",
                                            control = list(niter=10)))

    phy.unit <- as.unit.tree(fit.k)

    ## Manual rescaling using Geiger functions
    cmp <- rescale(phy, "kappa", kappa=coef(fit.k)[[1]], sigsq=coef(fit.k)[[2]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})

test_that("delta tree rescaling worked (fitContinuous)", {
    fit.d <- suppressWarnings(fitContinuous(phy=phy, dat=states, model="delta",
                                            control = list(niter=10)))

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
    fit.w <- suppressWarnings(fitContinuous(phy=phy, dat=states, model="white",
                                            control = list(niter=10)))

    phy.unit <- as.unit.tree(fit.w)

    ## Manual rescaling using Geiger function
    cmp <- rescale(phy, "white", sigsq=coef(fit.w)[[1]])

    expect_that(phy.unit$phy, equals(cmp))
    expect_that(is.unit.tree(phy.unit), is_true())
})
                
    
