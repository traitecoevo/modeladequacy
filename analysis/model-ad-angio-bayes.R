## Script for evaluating model adequacy of angiosperm clades
## Fitting models using maximum likelihood

## Read in all the helper files
source("R/modelfit-helper-fxns.R") 

## Load in data
all.dat <- get.angio.data()



## Wrapper function
## Fits BM, OU and EB using MCMC
## Computes DIC scores for all the models
## Assess model adequacy
##
## Function takes one of the elements of angio-trait-data-all.rds
model.ad.bayes <- function(x){

    ## Define objects
    phy    <- x$phy
    states <- x$states
    SE     <- x$SE
    idx    <- x$index
    taxa   <- x$taxa
    trait  <- x$trait
    rank   <- x$rank
    n      <- Ntip(phy)
    t      <- max(branching.times(phy))
    
    
    ## Create likelihood fxns
    lik.bm <- make.bm(phy, states, SE, control=dt.con)
    lik.ou <- make.ou(phy, states, SE, control=dt.con)
    lik.eb <- make.eb(phy, states, SE, control=dt.con)

    ## Start at REML estimate of sigsq for all models
    ## use function sigst.est to estimate this
    s2 <- sigsq.est(as.unit.tree(phy, data=states))

    ## define priors for analysis
    ## Prior for BM model
    prior.bm <- make.prior.bm(s2.lower=0, s2.upper=2)

    ## Prior for OU model
    prior.ou <- make.prior.ou(s2.lower=0, s2.upper=2,
                          ln.mean=log(0.5), ln.sd=log(1.5))

    ## Prior for EB model
    prior.eb <- make.prior.eb(s2.lower=0, s2.upper=2,
                          a.lower=-1, a.upper=0)


    ## Run short chain to obtain appropriate step size for MCMC
    tmp.bm <- mcmc(lik.bm, x.init=s2, nsteps=100,
                   prior=prior.bm, w=1, print.every=0)
    tmp.ou <- mcmc(lik.ou, x.init=c(s2, 0.05), nsteps=100,
                   prior=prior.ou, w=1, print.every=0)
    tmp.eb <- mcmc(lik.eb, x.init=c(s2, -0.1), nsteps=100,
                   prior=prior.eb, w=1, print.every=0)

    ## Calculate reasonable step size
    w.bm <- diff(sapply(tmp.bm[2], range))
    w.ou <- diff(sapply(tmp.ou[2:3], range))
    w.eb <- diff(sapply(tmp.eb[2:3], range))

    ## run MCMC chain for 10000 steps
    samp.bm <- mcmc(lik.bm, x.init=s2, nsteps=10000,
                    prior=prior.bm, w=w.bm, print.every=0)
    samp.ou <- mcmc(lik.ou, x.init=c(s2,0.05), nsteps=10000,
                    prior=prior.ou, w=w.ou, print.every=0)
    samp.eb <- mcmc(lik.eb, x.init=c(s2,-0.1), nsteps=10000,
                    prior=prior.eb, w=w.eb, print.every=0)

    ## Calculate dic and dic weights
    ## Remove 1000 samples as burnin
    dic.bm <- dic.mcmcsamples(samp.bm, burnin=1000)
    dic.ou <- dic.mcmcsamples(samp.ou, burnin=1000)
    dic.eb <- dic.mcmcsamples(samp.eb, burnin=1000)

    dic.all <- c(dic.bm, dic.ou, dic.eb)
    names(dic.all) <- c("BM", "OU", "EB")

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


    
    ## Create bayesian output
    out.bayes <- c(x$taxa, x$rank, x$trait, n, t,
             dic.bm, unname(dic.w["BM"]), pv.mcmc.bm, mv.mcmc.bm,
             dic.ou, unname(dic.w["OU"]), pv.mcmc.ou, mv.mcmc.ou,
             dic.eb, unname(dic.w["EB"]), pv.mcmc.eb, mv.mcmc.eb)

    out.bayes <- as.data.frame(out.bayes)
    rownames(out.bayes) <- c("taxa", "rank", "trait", "size", "age",
                    "dic.bm", "dicw.bm", paste(ss, "mcmc.bm", sep="."), "mv.mcmc.bm",
                    "dic.ou", "dicw.ou", paste(ss, "mcmc.ou", sep="."), "mv.mcmc.ou",
                    "dic.eb", "dicw.eb", paste(ss, "mcmc.eb", sep="."), "mv.mcmc.eb")

    write.csv(out.bayes, file=paste(paste("output/allres", "bayes", trait, idx, sep="_"),
                             ".csv", sep=""))
    
}



## IMPORTANT NOTE:
## Occasionally, due to apparent problems in optimization,
## a dataset may fail or get stuck.
## For this reason, I created the output for each dataset as its own file
## rather than simply outputing the results and combining all the
## results with do.call() after the fact.
## This is not ideal or elegant -- it creates a lot files in the output folder
## but was a way of ensuring the analyses all worked.



## define number of cores for parallelization
max.cores <- 32


## Run analyses across all clades
all.res <- mclapply(all.dat, function(x) model.ad.bayes(x),
                    mc.cores=max.cores, mc.preschedule=FALSE)




## Process results

## Gather all the files together and produce a single output file
tmp <- dir("output")

## pull out only the file names with allres_ml
bay <- tmp[grep("allres_bayes", tmp)]

## read in all results individually and bind them together
bay.res <- lapply(bay, function(x) build.adequacy.results(x))
bay.res <- do.call(rbind, bay.res)

## make sure numeric results are not factors
bay.res[,c(4:ncol(bay.res))] <- sapply(c(4:ncol(bay.res)), function(x)
                                     unfactor(bay.res[,x]))

## write to csv
write.csv(bay.res, file="output/bayes-results.csv")


