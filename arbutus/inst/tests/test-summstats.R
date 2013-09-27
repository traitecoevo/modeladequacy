source("helper-arbutus.R")

context("Summary statistics")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[,"wingL"]
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



test_that("REML sigsq is being calculated correctly", {
    phy.unit <- as.unit.tree(phy, states)

    ss <- round(sigsq.reml(phy.unit), 7)
    expect_that(ss, equals(0.130356))
})



    
    

                

    
