## # Analysis of model adequacy for angiosperm functional traits

## Get functions for reading in data
source("R/modelfit-helper-fxns.R") 

## # Pull in all of the information for all clades and traits

## Load in data
all.dat <- get.angio.data()

## # Create functions for fitting models and examining adequacy of each model

## For the sake of clarity, define a list which can
## passed to the 'control' argument of all diversitree's fxns
dt.con <- list(method="pruning", backend="C")



## Define function for 

## Define prior functions for mcmc analyses:
## BM prior
prior.bm <- make.prior.bm(s2.lower=0, s2.upper=2)

## OU prior
prior.ou <- make.prior.ou(s2.lower=0, s2.upper=2,
                          ln.mean=log(0.5), ln.sd=log(1.5))

## EB prior
prior.eb <- make.prior.eb(s2.lower=0, s2.upper=2,
                          a.lower=-1, a.upper=0)

model.ad.wrap <- function(x){

    ## Define objects
    phy <- x$phy
    states <- x$states
    SE <- x$SE
    
    ## Create likelihood fxns
    lik.bm <- make.bm(phy, states, SE, control=dt.con)
    lik.ou <- make.ou(phy, states, SE, control=dt.con)
    lik.eb <- make.eb(phy, states, SE, control=dt.con)


    ## Maximum likelihood
    ## Fit models
    fit.bm <- find.mle(lik.bm, x.init=0.1)
    fit.ou <- find.mle(lik.ou, x.init=c(0.1, 0.1))
    fit.eb <- find.mle(lik.eb, x.init=c(0.1, -0.1))

    ## Calculate aic and aic weights for all models
    aic.bm <- AIC(fit.bm)
    aic.ou <- AIC(fit.ou)
    aic.eb <- AIC(fit.eb)

    aic.all <- c(aic.bm, aic.ou, aic.eb)
    names(aic.all) <- c("BM", "OU", "EB")

    aic.w <- ic.weights(aic.all)

    ## Assess adequacy of all models
    ml.bm <- phy.model.check(fit.bm)
    ml.ou <- phy.model.check(fit.ou)
    ml.eb <- phy.model.check(fit.eb)

    ## Get p values for all summary stats
    pv.ml.bm <- as.numeric(pval.summ.stats(ml.bm))
    pv.ml.ou <- as.numeric(pval.summ.stats(ml.ou))
    pv.ml.eb <- as.numeric(pval.summ.stats(ml.eb))

    ## Get malahanobis distance for all summary stats
    mv.ml.bm <- mv.summ.stats(ml.bm)
    mv.ml.ou <- mv.summ.stats(ml.ou)
    mv.ml.eb <- mv.summ.stats(ml.eb)



    ## Bayesian inference
    ## Run short chain to obtain appropriate step size for MCMC
    ## Start with ML estimates
    tmp.bm <- mcmc(lik.bm, x.init=fit.bm$par, nsteps=100,
                   prior=prior.bm, w=1, print.every=0)
    tmp.ou <- mcmc(lik.ou, x.init=fit.ou$par, nsteps=100,
                   prior=prior.ou, w=1, print.every=0)
    tmp.eb <- mcmc(lik.eb, x.init=fit.eb$par, nsteps=100,
                   prior=prior.eb, w=1, print.every=0)

    ## Calculate reasonable step size
    w.bm <- diff(sapply(tmp.bm[2], range))
    w.ou <- diff(sapply(tmp.ou[2:3], range))
    w.eb <- diff(sapply(tmp.eb[2:3], range))

    ## run MCMC chain for 10000 steps
    samp.bm <- mcmc(lik.bm, x.init=fit.bm$par, nsteps=10000,
                    prior=prior.bm, w=w.bm, print.every=0)
    samp.ou <- mcmc(lik.ou, x.init=fit.ou$par, nsteps=10000,
                    prior=prior.ou, w=w.ou, print.every=0)
    samp.eb <- mcmc(lik.eb, x.init=fit.eb$par, nsteps=10000,
                    prior=prior.eb, w=w.eb, print.every=0)

    ## Calculate dic and dic weights
    ## Remove 1000 samples as burnin
    dic.bm <- dic.mcmcsamples(samp.bm, burnin=1000)
    dic.ou <- dic.mcmcsamples(samp.ou, burnin=1000)
    dic.eb <- dic.mcmcsamples(samp.eb, burnin=1000)

    dic.all <- c(dic.bm, dic.ou, dic.eb)
    names(dic.all) <- names(aic.all)

    dic.w <- ic.weights(dic.all)
    
    ## Assess adequacy of all models
    mcmc.bm <- phy.model.check(samp.bm, burnin=1000, sample=1000)
    mcmc.ou <- phy.model.check(samp.ou, burnin=1000, sample=1000)
    mcmc.eb <- phy.model.check(samp.eb, burnin=1000, sample=1000)

    ## Get p values for all summary stats
    pv.mcmc.bm <- as.numeric(pval.summ.stats(mcmc.bm))
    pv.mcmc.ou <- as.numeric(pval.summ.stats(mcmc.ou))
    pv.mcmc.eb <- as.numeric(pval.summ.stats(mcmc.eb))

    ## Get malahanobis distance for all summary stats
    mv.mcmc.bm <- mean(mv.summ.stats(mcmc.bm))
    mv.mcmc.ou <- mean(mv.summ.stats(mcmc.ou))
    mv.mcmc.eb <- mean(mv.summ.stats(mcmc.eb))


    
    ## Create output
    n <- Ntip(phy)
    t <- max(branching.times(phy))
    out <- c(x$taxa, x$rank, x$trait, n, t,
             aic.bm, unname(aic.w["BM"]), pv.ml.bm, mv.ml.bm,
             aic.ou, unname(aic.w["OU"]), pv.ml.ou, mv.ml.ou,
             aic.eb, unname(aic.w["EB"]), pv.ml.eb, mv.ml.eb,
             dic.bm, unname(dic.w["BM"]), pv.mcmc.bm, mv.mcmc.bm,
             dic.ou, unname(dic.w["OU"]), pv.mcmc.ou, mv.mcmc.ou,
             dic.eb, unname(dic.w["EB"]), pv.mcmc.eb, mv.mcmc.eb)

    ss <- names(def.summ.stats())
    names(out) <- c("taxa", "rank", "trait", "size", "age",
                    "aic.bm", "aicw.bm", paste(ss, "ml.bm", sep="."), "mv.ml.bm",
                    "aic.ou", "aicw.ou", paste(ss, "ml.ou", sep="."), "mv.ml.ou",
                    "aic.eb", "aicw.eb", paste(ss, "ml.eb", sep="."), "mv.ml.eb",
                    "dic.bm", "dicw.bm", paste(ss, "mcmc.bm", sep="."), "mv.mcmc.bm",
                    "dic.ou", "dicw.ou", paste(ss, "mcmc.ou", sep="."), "mv.mcmc.ou",
                    "dic.eb", "dicw.eb", paste(ss, "mcmc.eb", sep="."), "mv.mcmc.eb")

    ## Return output
    out
}







## # Run analyses across all clades

## define number of cores for parallelization
max.cores <- 3
all.res <- mclapply(all.dat, function(x) model.ad.wrap(x),
                    mc.cores=max.cores, mc.preschedule=FALSE)

res <- do.call(rbind, all.res)

## write csv with all resutls
write.csv(res, "output/all-res-angio.csv")





## # Plotting functions HERE 



