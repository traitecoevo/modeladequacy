source("helper-arbutus.R")

context("Tree rescaling functions")

se <- 0.01
data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- dat$phy

f <- function(p) c(as.list(p), SE=0)
g <- function(p) c(as.list(p), SE=se)

test_that("Rescaling BM agrees with geiger", {
  p0.bm <- c(sigsq=1)    # no change
  p.bm <- c(sigsq=0.06) # rescale

  expect_that(arbutus:::model.phylo.bm(phy, f(p0.bm)),
              is_identical_to(phy))
  expect_that(arbutus:::model.phylo.bm(phy, f(p.bm)),
              equals(rescale(phy, "BM", p.bm)))
  expect_that(arbutus:::model.phylo.bm(phy, g(p.bm)),
              equals(rescale.se(phy, "BM", p.bm, SE=se)))
})

test_that("Rescaling OU agrees with geiger", {
  p0.ou <- c(alpha=0, sigsq=1)
  p1.ou <- c(alpha=2.7, sigsq=0.11)
  p2.ou <- c(alpha=1e-20, sigsq=1)

  expect_that(arbutus:::model.phylo.ou(phy, f(p0.ou)),
              is_identical_to(phy))

  expect_that(arbutus:::model.phylo.ou(phy, f(p1.ou)),
              equals(rescale(phy, "OU", p1.ou[[1]], p1.ou[[2]])))
  expect_that(arbutus:::model.phylo.ou(phy, f(p2.ou)),
              equals(phy))

  expect_that(arbutus:::model.phylo.ou(phy, g(p1.ou)),
              equals(rescale.se(phy, "OU", p1.ou[[1]], p1.ou[[2]], SE=se)))
  expect_that(arbutus:::model.phylo.ou(phy, g(p2.ou)),
              equals(rescale.se(phy, "BM", p2.ou["sigsq"], SE=se)))
})

test_that("Rescaling EB agrees with geiger", {
  p0.eb <- c(a=0, sigsq=1)
  p.eb <- c(a=-0.01, sigsq=0.06)

  expect_that(arbutus:::model.phylo.eb(phy, f(p0.eb)),
              is_identical_to(phy))
  expect_that(arbutus:::model.phylo.eb(phy, f(p.eb)),
              equals(rescale(phy, "EB", p.eb[[1]], p.eb[[2]])))
  expect_that(arbutus:::model.phylo.eb(phy, g(p.eb)),
              equals(rescale.se(phy, "EB", p.eb[[1]], p.eb[[2]], SE=se)))
})

test_that("Rescaling lambda agrees with geiger", {
  p0.lambda <- c(lambda=1, sigsq=1)
  p.lambda <- c(lambda=.7, sigsq=0.01)

  expect_that(arbutus:::model.phylo.lambda(phy, f(p0.lambda)),
              is_identical_to(phy))
  expect_that(arbutus:::model.phylo.lambda(phy, f(p.lambda)),
              equals(rescale(phy, "lambda", p.lambda[[1]], p.lambda[[2]])))
  expect_that(arbutus:::model.phylo.lambda(phy, g(p.lambda)),
              equals(rescale.se(phy, "lambda", p.lambda[[1]],
                                p.lambda[[2]], SE=se)))
})

test_that("Rescaling kappa aggrees with geiger", {
  p0.kappa <- c(kappa=1, sigsq=1)
  p.kappa <- c(kappa=0.6, sigsq=0.03)

  expect_that(arbutus:::model.phylo.kappa(phy, f(p0.kappa)),
              is_identical_to(phy))
  expect_that(arbutus:::model.phylo.kappa(phy, f(p.kappa)),
              equals(rescale(phy, "kappa", p.kappa[[1]], p.kappa[[2]])))
  expect_that(arbutus:::model.phylo.kappa(phy, g(p.kappa)),
              equals(rescale.se(phy, "kappa", p.kappa[[1]], p.kappa[[2]],
                                SE=se)))
})

## TODO: Matt, please look into this:
test_that("Rescaling delta agrees with geiger", {
  p0.delta <- c(delta=1, sigsq=1)
  p.delta <- c(delta=3, sigsq=0.04)

  expect_that(arbutus:::model.phylo.delta(phy, f(p0.delta)),
              equals(phy)) # does transform
  ## expect_that(arbutus:::model.phylo.delta(phy, f(p.delta)),
  ##             equals(rescale(phy, "delta", p.delta[[1]], p.delta[[2]])))
  ## expect_that(arbutus:::model.phylo.delta(phy, g(p.delta)),
  ##             equals(rescale.se(phy, "delta", p.delta[[1]], p.delta[[2]],
  ##                            SE=se)))
})

test_that("Rescaling white noise agrees with geiger", {
  p.white <- c(sigsq=0.03) # there is no non-transform case here.
  expect_that(arbutus:::model.phylo.white(phy, f(p.white)),
              equals(rescale(phy, "white", p.white)))
  expect_that(arbutus:::model.phylo.white(phy, g(p.white)),
              equals(rescale.se(phy, "white", p.white, SE=se)))
})
