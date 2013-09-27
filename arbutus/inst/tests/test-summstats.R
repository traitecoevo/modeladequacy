source("helper-arbutus.R")

context("Summary statistics")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[,"wingL"]

test_that("Summary statistic wrapper function produces correct output", {
    phy.unit <- as.unit.tree(phy, states)

    ## check default summary stats
    stat.fxn <- arbutus:::defSummStats()
    expect_that(is.list(stat.fxn), is_true())
    expect_that(length(stat.fxn), equals(6))

    ## make sure it uses correct summary stats if default
    ss <- summStat(phy.unit, stats=NULL)
    ss.names <- names(ss)
    expect_that(ss.names, equals(names(stat.fxn)))

    ## if custom fxn supplied
    foo <- function(x)
        mean(x$pics[,"contrasts"])

    ss.mean <- summStat(phy.unit, stats=list(mean=foo))
    ss.names.mean <- names(ss.mean)
    expect_that(ss.names.mean, equals("mean"))

    ## make sure output is correct for a single and multiple unit trees
    expect_that(ss, is_a("data.frame"))
    nrow.ss <- nrow(ss)
    expect_that(nrow.ss, equals(1))
    ut <- list(phy.unit, phy.unit)
    set.seed(1)
    ss.multi <- summStat(ut, stats=NULL)
    expect_that(ss.multi, is_a("data.frame"))
    nrow.ss.multi <- nrow(ss.multi)
    expect_that(nrow.ss.multi, equals(2))

})

test_that("REML sigsq is being calculated correctly", {
    phy.unit <- as.unit.tree(phy, states)

    ss <- round(as.numeric(summStat(phy.unit)[,"reml.sigsq"]), 7)
    expect_that(ss, equals(0.130356))
})



    
    

                

    
