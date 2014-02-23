source("helper-arbutus.R")

context("Summary statistics")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[phy$tip.label,"wingL"]
phy.unit <- as.unit.tree(phy, states)

test_that("Summary statistic wrapper function produces correct output", {
    ## check default summary stats
    stat.fxn <- arbutus:::def.summ.stats()
    expect_that(is.list(stat.fxn), is_true())
    expect_that(length(stat.fxn), equals(6))

    ## make sure it uses correct summary stats if default
    ss <- summ.stats(phy.unit, stats=NULL)
    ss.names <- names(ss)
    expect_that(ss.names, equals(names(stat.fxn)))

    ## Single tree:
    expect_that(ss, is_a("data.frame"))
    expect_that(nrow(ss), equals(1))

    ## Multiple trees:
    ut <- list(phy.unit, phy.unit)
    set.seed(1)
    ss.multi <- summ.stats(ut, stats=NULL)
    expect_that(ss.multi, is_a("data.frame"))
    nrow.ss.multi <- nrow(ss.multi)
    expect_that(nrow.ss.multi, equals(2))

    ## Something else
    expect_that(summ.stats(list()), throws_error())
    expect_that(summ.stats(phy), throws_error())
    expect_that(summ.stats(1), throws_error())
})

test_that("Custom supplied functions work correctly", {
    ## if custom fxn supplied
    foo <- function(x)
        mean(x$pics[,"contrasts"])

    ss.mean <- summ.stats(phy.unit, stats=list(mean=foo))
    ss.names.mean <- names(ss.mean)
    expect_that(ss.names.mean, equals("mean"))

    ## Things that should cause errors:
    ## Unnamed list:
    expect_that(summ.stats(phy.unit, stats=list(foo)),
                throws_error())
    ## Not given a list
    expect_that(summ.stats(phy.unit, stats=foo),
                throws_error())
    ## Empty list
    expect_that(summ.stats(phy.unit, stats=list()),
                throws_error())
    ## Not functions
    expect_that(summ.stats(phy.unit, stats=list(a=1)),
                throws_error())
})



## Numbers here computed using make_expectation()
test_that("REML sigsq is being calculated correctly", {
  expect_that(sigsq.est(phy.unit), equals(0.0764244853752198))
})

test_that("Coefficient of variation of contrasts is being calculated correctly", {
  expect_that(cv.contrast(phy.unit), equals(0.7362652135189863))
})

test_that("Slope of contrasts and variance is being calculated correctly", {
  expect_that(cor.contrast.var(phy.unit), equals(-0.29076711665925))
})

test_that("Slope of contrasts and asr is being calculated correctly", {
  expect_that(cor.contrast.asr(phy.unit), equals(-1.22336966244042))
})

test_that("Slope of contrasts and node height is being calculated correctly", {
  expect_that(cor.contrast.nh(phy.unit), equals(-0.306687992622041))
})

test_that("KS D-Statistic is being calculated correctly",{
    set.seed(1) ## simulates normal so need to set seed
    expect_that(ks.contrast(phy.unit), equals(0.198466666666667))
})


test_that("Summary statistic comparison is working properly",{
    ## single observed summary stat
    set.seed(1)
    ss.o <- summ.stats(phy.unit)
    stats <- names(ss.o)

    set.seed(1)
    sims <- sim.char.unit(phy.unit, nsim=6)

    set.seed(1)
    ss.s <- summ.stats(sims)
    cc <- compare.summ.stats(ss.o, ss.s)

    expect_that(cc, is_a("phy.ss"))
    expect_that(cc, is_a("list"))
    el <- c("p.values", "summ.stats.obs", "summ.stats.sim")
    expect_that(names(cc), equals(el))
    expect_that(cc$p.values, is_a("numeric"))
    expect_that(cc$summ.stats.obs, is_a("data.frame"))
    expect_that(cc$summ.stats.sim, is_a("data.frame"))
    expect_that(names(cc$p.values), equals(stats))
    expect_that(names(cc$summ.stats.obs), equals(stats))
    expect_that(names(cc$summ.stats.sim), equals(stats))
    expect_that(nrow(cc$summ.stats.obs), equals(1))
    expect_that(nrow(cc$summ.stats.sim), equals(6))

    ## check if values are correct
    expect_that(as.numeric(cc$p.values),
                equals(c(0, 0, 0.285714285714286,
                         0, 0.285714285714286, 0.285714285714286)))

    ## if we provide different names
    ss.s.wrgname <- ss.s
    names(ss.s.wrgname) <- c("a", "b", "c", "d", "e", "f")
    expect_that(compare.summ.stats(ss.o, ss.wrgname), throws_error())

    ## if lengths of summary stats are not the same
    ss.wrgnum <- ss.s[c(1:5)]
    expect_that(compare.summ.stats(ss.o, ss.wrgnum), throws_error())


    ## check for multiple datasets
    set.seed(1)
    ss.o <- cbind.data.frame(rnorm(6), rnorm(6), rnorm(6),
                             rnorm(6), rnorm(6), rnorm(6))
    set.seed(42)
    ss.s <- cbind.data.frame(rnorm(6), rnorm(6), rnorm(6),
                             rnorm(6), rnorm(6), rnorm(6))
    names(ss.o) <- names(ss.s) <- stats

    cm <- compare.summ.stats(ss.o, ss.s)
    expect_that(cm, is_a("phy.ss"))
    expect_that(cm, is_a("list"))
    el <- c("p.values", "summ.stats.obs", "summ.stats.sim")
    expect_that(names(cm), equals(el))
    expect_that(cm$p.values, is_a("numeric"))
    expect_that(cm$summ.stats.obs, is_a("data.frame"))
    expect_that(cm$summ.stats.sim, is_a("data.frame"))
    expect_that(names(cm$p.values), equals(stats))
    expect_that(names(cm$summ.stats.obs), equals(stats))
    expect_that(names(cm$summ.stats.sim), equals(stats))
    expect_that(nrow(cm$summ.stats.obs), equals(6))
    expect_that(nrow(cm$summ.stats.sim), equals(6))

    ex.p <- c(4/6, 4/6, 4/6, 4/6, 4/6, 1)
    expect_that(as.numeric(cm$p.values), equals(ex.p))

    ## if we provide different names
    ss.s.wrgname <- ss.s
    names(ss.s.wrgname) <- c("a", "b", "c", "d", "e", "f")
    expect_that(compare.summ.stats(ss.o, ss.wrgname), throws_error())

    ## if lengths of summary stats are not the same
    ss.wrgnum <- ss.s[,c(1:5)]
    expect_that(compare.summ.stats(ss.o, ss.wrgnum), throws_error())

    ## if number of rows are not the same
    ss.wrgrow <- ss.s[c(1:5),]
    expect_that(compare.summ.stats(ss.o, ss.wrgnum), throws_error())

})


test_that("p values are properly extracted from phy.ss object",{
    ss.o <- summ.stats(phy.unit)
    stats <- names(ss.o)
    sims <- sim.char.unit(phy.unit, nsim=6)
    ss.s <- summ.stats(sims)
    cc <- compare.summ.stats(ss.o, ss.s)

    pv <- cc$p.values
    expect_that(pval.summ.stats(cc), equals(pv))
})

test_that("mahalanobis distance is calculated properly from summary stats",{
    ## single set of summary stats supplied
    ss.o <- summ.stats(phy.unit)
    ## need to 100 to have enough info for covariance
    sims <- sim.char.unit(phy.unit, nsim=100)
    ss.s <- summ.stats(sims)
    cc <- compare.summ.stats(ss.o, ss.s)

    mv <- mv.summ.stats(cc)
    expect_that(length(mv), equals(1))

    ## calculated by hand
    ## log d.ks
    cc$summ.stats.obs[,"d.ks"] <- log(cc$summ.stats.obs[,"d.ks"])
    cc$summ.stats.sim[,"d.ks"] <- log(cc$summ.stats.sim[,"d.ks"])

    obs.ss <- as.matrix(cc$summ.stats.obs)
    sim.ss <- as.matrix(cc$summ.stats.sim)
    cv <- cov(sim.ss)

    m <- mahalanobis(x=obs.ss, center=colMeans(sim.ss), cov=cv)

    expect_that(mv, equals(m))

    ## multiple set of summary statitics supplied
    m.unit <- lapply(c(1:100), function(x) return(phy.unit))
    ss.o <- summ.stats(m.unit)
    sims <- sim.char.unit(m.unit)
    ss.s <- summ.stats(sims)
    cc <- compare.summ.stats(ss.o, ss.s)

    mv <- mv.summ.stats(cc)
    expect_that(length(mv), equals(100))

    ## calculated by hand
    ## log ks.contrast
    cc$summ.stats.obs[,"d.ks"] <- log(cc$summ.stats.obs[,"d.ks"])
    cc$summ.stats.sim[,"d.ks"] <- log(cc$summ.stats.sim[,"d.ks"])

    obs.ss <- as.matrix(cc$summ.stats.obs)
    sim.ss <- as.matrix(cc$summ.stats.sim)
    cv <- cov(sim.ss)

    m <- mahalanobis(x=obs.ss, center=colMeans(sim.ss), cov=cv)

    expect_that(mv, equals(m))

    ## check to make sure that if a non phy.ss object supplied
    ## mv.summ.stats breaks
    expect_that(mv.summ.stats(list()), throws_error())
})

    
    

                

    
