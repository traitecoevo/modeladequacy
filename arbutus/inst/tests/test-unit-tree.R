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
