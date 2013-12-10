source("helper-arbutus.R")

context("Tree rescaling functions")

test_that("Rescaling functions work", {
  p.bm <- c(sigsq=0.06)
  p.ou <- c(alpha=2.7, sigsq=0.11)
  p.eb <- c(a=-0.01, sigsq=0.06)
  p.la <- c(lambda=0.5, sigsq=0.01)
  p.ka <- c(kappa=0.6, sigsq=0.03)
  p.de <- c(delta=3, sigsq=0.04)
  p.wh <- c(sigsq=0.03)
  se <- 0.01

  data(geospiza)
  dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

  phy <- dat$phy
  states <- dat$dat[,"wingL"]

  f <- function(p) c(as.list(p), SE=0)
  g <- function(p) c(as.list(p), SE=se)

  expect_that(arbutus:::model.phylo.bm(phy, f(p.bm)),
              equals(rescale(phy, "BM", p.bm)))
  expect_that(arbutus:::model.phylo.ou(phy, f(p.ou)),
              equals(rescale(phy, "OU", p.ou[[1]], p.ou[[2]])))
  expect_that(arbutus:::model.phylo.eb(phy, f(p.eb)),
              equals(rescale(phy, "EB", p.eb[[1]], p.eb[[2]])))
  expect_that(arbutus:::model.phylo.lambda(phy, f(p.la)),
              equals(rescale(phy, "lambda", p.la[[1]], p.la[[2]])))
  expect_that(arbutus:::model.phylo.kappa(phy, f(p.ka)),
              equals(rescale(phy, "kappa", p.ka[[1]], p.ka[[2]])))
  ## TODO: Matt - please look here
  ## expect_that(arbutus:::model.phylo.delta(phy, f(p.de)),
  ##             equals(rescale(phy, "delta", p.de[[1]], p.de[[2]])))
  expect_that(arbutus:::model.phylo.white(phy, f(p.wh)),
              equals(rescale(phy, "white", p.wh)))

  expect_that(arbutus:::model.phylo.bm(phy, g(p.bm)),
              equals(rescale.se(phy, "BM", p.bm, SE=se)))
  expect_that(arbutus:::model.phylo.ou(phy, g(p.ou)),
              equals(rescale.se(phy, "OU", p.ou[[1]], p.ou[[2]], SE=se)))
  expect_that(arbutus:::model.phylo.eb(phy, g(p.eb)),
              equals(rescale.se(phy, "EB", p.eb[[1]], p.eb[[2]], SE=se)))
  expect_that(arbutus:::model.phylo.lambda(phy, g(p.la)),
              equals(rescale.se(phy, "lambda", p.la[[1]], p.la[[2]], SE=se)))
  expect_that(arbutus:::model.phylo.kappa(phy, g(p.ka)),
              equals(rescale.se(phy, "kappa", p.ka[[1]], p.ka[[2]], SE=se)))
  ## TODO: Matt - please look here
  ## expect_that(arbutus:::model.phylo.delta(phy, g(p.de)),
  ##             equals(rescale.se(phy, "delta", p.de[[1]], p.de[[2]], SE=se)))
  expect_that(arbutus:::model.phylo.white(phy, g(p.wh)),
              equals(rescale.se(phy, "white", p.wh, SE=se)))
})
