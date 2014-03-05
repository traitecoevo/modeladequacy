## Script for evaluating model adequacy of angiosperm clades
## Fitting models using maximum likelihood

## Read in all the helper files
source("R/modelfit-helper-fxns.R") 

## Load in data
all.dat <- get.angio.data()


## Wrapper function
## Fits BM, OU and EB
## Computes AIC scores for all the models
## Assess model adequacy
##
## Function takes one of the elements of angio-trait-data-all.rds
model.ad.ml <- function(x){

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
    s2 <- sigst.est(as.unit.tree(phy, data=states))
    
    ## get upper bounds for OU model
    upbnd.ou <- bnds.ou(phy, states)
    
    ## Fit models
    fit.bm <- find.mle(lik.bm, x.init=s2)
    fit.ou <- find.mle(lik.ou, x.init=c(s2, 0.05), upper=upbnd.ou)
    fit.eb <- find.mle(lik.eb, x.init=c(s2, -0.1), lower=c(0,-1))

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

    ## Mean diagnostic
    ## This is not strictly necessary but was used to
    ## test whether m.pic was consistently under or overestimated
    ## This is the only summary statistic which we expected to be
    ## systematically biased in one direction
    md.bm <- length(which(ml.bm$summ.stats.sim[,"m.pic"] < ml.bm$summ.stats.obs[,"m.pic"]))
    md.ou <- length(which(ml.ou$summ.stats.sim[,"m.pic"] < ml.ou$summ.stats.obs[,"m.pic"]))
    md.eb <- length(which(ml.eb$summ.stats.sim[,"m.pic"] < ml.eb$summ.stats.obs[,"m.pic"]))

    ## write out summary file
    out.ml <- c(x$taxa, x$rank, x$trait, n, t,
                aic.bm, unname(aic.w["BM"]), pv.ml.bm, mv.ml.bm, md.bm,
                aic.ou, unname(aic.w["OU"]), pv.ml.ou, mv.ml.ou, md.ou,
                aic.eb, unname(aic.w["EB"]), pv.ml.eb, mv.ml.eb, md.eb)

    out.ml <- as.data.frame(out.ml)
    rownames(out.ml) <- c("taxa", "rank", "trait", "size", "age",
                       "aic.bm", "aicw.bm", paste(ss, "ml.bm", sep="."), "mv.ml.bm", "mean.diag.bm",
                       "aic.ou", "aicw.ou", paste(ss, "ml.ou", sep="."), "mv.ml.ou", "mean.diag.ou",
                       "aic.eb", "aicw.eb", paste(ss, "ml.eb", sep="."), "mv.ml.eb", "mean.diag.eb")

    write.csv(out.ml, file=paste(paste("output/allres", "ml", trait, idx, sep="_"), ".csv", sep=""))
}


## IMPORTANT NOTE:
## Occasionally, due to apparent problems in optimization,
## a dataset may fail or get stuck.
## For this reason, I created the output for each dataset as its own file
## rather than simply outputing the results and combining all the
## results with do.call() after the fact.
## This is not ideal or elegant -- it creates a lot files in the output folder
## but was a way of ensuring the analyses all worked.
##
## When the analyses were run for the paper, there were datasets that
## had problems. I went back and ran each of these "by hand", by altering starting
## parameters and optimizer conditions in order to ensure convergence.




## define number of cores for parallelization
max.cores <- 32


## Run analyses across all clades
all.res <- mclapply(all.dat, function(x) model.ad.ml(x),
                    mc.cores=max.cores, mc.preschedule=FALSE)




## Process results

## Gather all the files together and produce a single output file
tmp <- dir("output")

## pull out only the file names with allres_ml
ml <- tmp[grep("allres_ml", tmp)]

## function for reading in individual ml results
build.adequacy.results.ml <- function(x){
    f <- read.csv(paste("output", x, sep="/"), as.is=TRUE, row.names=1)
    d <- as.data.frame(t(f))
    rownames(d) <- NULL
    d
}

## read in all results individually and bind them together
ml.res <- lapply(ml, function(x) build.adequacy.results.ml(x))
ml.res <- do.call(rbind, ml.res)

## unfactor all numeric columns
unfactor <- function(x)
    as.numeric(levels(x))[x]

ml.res[,c(4:ncol(ml.res))] <- sapply(c(4:ncol(ml.res)), function(x)
                                     unfactor(ml.res[,x]))

## write to csv
write.csv(ml.res, file="output/ml-results.csv")



