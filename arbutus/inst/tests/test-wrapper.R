source("helper-arbutus.R")

context("wrapper (phy.model.check)")

data(finch)
phy <- finch$phy
dat1 <- finch$data[,"wingL"]
data <- as.data.frame(finch$data)
multiphy <- list(phy, phy)
class(multiphy) <- "multiPhylo"


## fit bm using fitContinuous
fitc.res <- fitContinuous(phy, dat=dat1, SE=0, control=list(niter=10))


## fit bm using diversitree
bm.lik <- make.bm(phy, dat1, control=list(method="pruning", backend="C"))

## mle
bm.mle.res <- find.mle(bm.lik, x.init = 1)

## mcmc
bm.mcmc.res <- mcmc(bm.lik, x.init = 1, nsteps=10, w=1, print.every=0)


## fit correlation with gls
gls.res <- gls(wingL~tarsusL, data=data, correlation = corBrownian(phy=phy))


## fit correlation with caper
dd <- cbind(data, phy$tip.label)
colnames(dd)[ncol(dd)] <- "species"
cd <- comparative.data(phy, dd, names.col = species)

caper.res <- pgls(wingL~tarsusL, data=cd)

## fit correlation with phylolm
phylolm.res <- phylolm(wingL~tarsusL, data=data, phy=phy, model="BM")


## fit correlation with pgls diversitree
pgls.lik <- make.pgls(phy, wingL~tarsusL, data=data)

## find mle for pgls model
pgls.mle.res <- find.mle(pgls.lik, x.init = c(0.5, 0.5, 0.5))

pgls.mcmc.res <- mcmc(pgls.lik, x.init = c(0.5, 0.5, 0.5), nsteps=10, w=1, print.every = 0)


test_that("phy model check runs properly", {
    mc.fitc <- phy.model.check(fitc.res, nsim=10)
    mc.bm.mle <- phy.model.check(bm.mle.res, nsim=10)
    mc.bm.mcmc <- phy.model.check(bm.mcmc.res, nsim=10)
    mc.gls <- phy.model.check(gls.res, nsim=10)
    mc.caper <- phy.model.check(caper.res, nsim=10)
    mc.phylolm <- phy.model.check(phylolm.res, nsim=10, phy=phy)
    mc.pgls.mle <- phy.model.check(pgls.mle.res, nsim=10)
    mc.pgls.mcmc <- phy.model.check(pgls.mcmc.res, nsim=10)
    mc.phy <- phy.model.check(phy, data=dat1, nsim=10)
    mc.multiphy <- phy.model.check(multiphy, data=dat1, nsim=10)
    
    expect_that(mc.fitc, is_a("phy.ss"))
    expect_that(mc.bm.mle, is_a("phy.ss"))
    expect_that(mc.bm.mcmc, is_a("phy.ss"))
    expect_that(mc.gls, is_a("phy.ss"))
    expect_that(mc.caper, is_a("phy.ss"))
    expect_that(mc.phylolm, is_a("phy.ss"))
    expect_that(mc.pgls.mle, is_a("phy.ss"))
    expect_that(mc.pgls.mcmc, is_a("phy.ss"))
    expect_that(mc.phy, is_a("phy.ss"))
    expect_that(mc.multiphy, is_a("phy.ss"))

    expect_that(nrow(mc.fitc$summ.stats.sim), equals(10))
    expect_that(nrow(mc.bm.mle$summ.stats.sim), equals(10))
    expect_that(nrow(mc.bm.mcmc$summ.stats.sim), equals(10))
    expect_that(nrow(mc.gls$summ.stats.sim), equals(10))
    expect_that(nrow(mc.caper$summ.stats.sim), equals(10))
    expect_that(nrow(mc.phylolm$summ.stats.sim), equals(10))
    expect_that(nrow(mc.pgls.mle$summ.stats.sim), equals(10))
    expect_that(nrow(mc.pgls.mcmc$summ.stats.sim), equals(10))
    expect_that(nrow(mc.phy$summ.stats.sim), equals(10))
    expect_that(nrow(mc.multiphy$summ.stats.sim), equals(2))

    expect_that(nrow(mc.fitc$summ.stats.obs), equals(1))
    expect_that(nrow(mc.bm.mle$summ.stats.obs), equals(1))
    expect_that(nrow(mc.bm.mcmc$summ.stats.obs), equals(10))
    expect_that(nrow(mc.gls$summ.stats.obs), equals(1))
    expect_that(nrow(mc.caper$summ.stats.obs), equals(1))
    expect_that(nrow(mc.phylolm$summ.stats.obs), equals(1))
    expect_that(nrow(mc.pgls.mle$summ.stats.obs), equals(1))
    expect_that(nrow(mc.pgls.mcmc$summ.stats.obs), equals(10))
    expect_that(nrow(mc.phy$summ.stats.obs), equals(1))
    expect_that(nrow(mc.multiphy$summ.stats.obs), equals(2))

    stats <- unlist(names(def.summ.stats()))

    st.fitc <- names(pval.summ.stats(mc.fitc))
    st.bm.mle <- names(pval.summ.stats(mc.bm.mle))
    st.bm.mcmc <- names(pval.summ.stats(mc.bm.mcmc))
    st.gls <- names(pval.summ.stats(mc.gls))
    st.caper <- names(pval.summ.stats(mc.caper))
    st.phylolm <- names(pval.summ.stats(mc.phylolm))
    st.pgls.mle <- names(pval.summ.stats(mc.pgls.mle))
    st.pgls.mcmc <- names(pval.summ.stats(mc.pgls.mcmc))
    st.phy <- names(pval.summ.stats(mc.phy))
    st.multiphy <- names(pval.summ.stats(mc.multiphy))

    expect_that(st.fitc, equals(stats))
    expect_that(st.bm.mle, equals(stats))
    expect_that(st.bm.mcmc, equals(stats))
    expect_that(st.gls, equals(stats))
    expect_that(st.caper, equals(stats))
    expect_that(st.phylolm, equals(stats))
    expect_that(st.pgls.mle, equals(stats))
    expect_that(st.pgls.mcmc, equals(stats))
    expect_that(st.phy, equals(stats))
    expect_that(st.multiphy, equals(stats))

})


test_that("if incomplete data is supplied, method fails",{
    
    ## phylolm without phylogeny
    expect_that(phy.model.check(phylolm.res), throws_error())

    ## phylogeny without data
    expect_that(phy.model.check(phy), throws_error())
    expect_that(phy.model.check(multiphy), throws_error())
})


test_that("stats are passed from wrapper to internals", {
    ## define some function
    foo <- function(x)
        mean(x$pics[,"contrasts"])
 

    ## create a named list with function
    new.stat <- list(mean=foo)
    
    mc.fitc <- phy.model.check(fitc.res, nsim=10, stats=new.stat)
    mc.bm.mle <- phy.model.check(bm.mle.res, nsim=10, stats=new.stat)
    mc.bm.mcmc <- phy.model.check(bm.mcmc.res, nsim=10, stats=new.stat)
    mc.gls <- phy.model.check(gls.res, nsim=10, stats=new.stat)
    mc.caper <- phy.model.check(caper.res, nsim=10, stats=new.stat)
    mc.phylolm <- phy.model.check(phylolm.res, nsim=10, stats=new.stat, phy=phy)
    mc.pgls.mle <- phy.model.check(pgls.mle.res, nsim=10, stats=new.stat)
    mc.pgls.mcmc <- phy.model.check(pgls.mcmc.res, nsim=10, stats=new.stat)
    mc.phy <- phy.model.check(phy, data=dat1, nsim=10, stats=new.stat)
    mc.multiphy <- phy.model.check(multiphy, data=dat1, nsim=10, stats=new.stat)

    st.fitc <- names(pval.summ.stats(mc.fitc))
    st.bm.mle <- names(pval.summ.stats(mc.bm.mle))
    st.bm.mcmc <- names(pval.summ.stats(mc.bm.mcmc))
    st.gls <- names(pval.summ.stats(mc.gls))
    st.caper <- names(pval.summ.stats(mc.caper))
    st.phylolm <- names(pval.summ.stats(mc.phylolm))
    st.pgls.mle <- names(pval.summ.stats(mc.pgls.mle))
    st.pgls.mcmc <- names(pval.summ.stats(mc.pgls.mcmc))
    st.phy <- names(pval.summ.stats(mc.phy))
    st.multiphy <- names(pval.summ.stats(mc.multiphy))

    stats <- "mean"
    
    expect_that(st.fitc, equals(stats))
    expect_that(st.bm.mle, equals(stats))
    expect_that(st.bm.mcmc, equals(stats))
    expect_that(st.gls, equals(stats))
    expect_that(st.caper, equals(stats))
    expect_that(st.phylolm, equals(stats))
    expect_that(st.pgls.mle, equals(stats))
    expect_that(st.pgls.mcmc, equals(stats))
    expect_that(st.phy, equals(stats))
    expect_that(st.multiphy, equals(stats))
})


    


    
    
