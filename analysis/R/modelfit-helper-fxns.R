## Functions used in analyses
## But not part of main analyses

## Load in all packages used 
library(arbutus)
library(diversitree)
library(multicore)

## Compute weights for AIC or DIC
## takes a named vector of AIC/DIC values
ic.weights <- function(x){
    
    d.x <- x - min(x)
    tmp <- exp(-0.5 * d.x)
    w   <- tmp / sum(tmp)
    names(w) <- names(x)
    w
}


## Compute deviance information criterion from mcmcsamples
dic.mcmcsamples <- function(x, burnin=NULL){
    if (!inherits(x, "mcmcsamples"))
        stop("this function is only designed for diversitrees mcmcsamples object")

    f <- attr(x, "func")

    if (!is.null(burnin)){
        if (burnin == 0)
            stop("if burnin supplied, it cannot be 0, leave as NULL")
        
        x <- x[-seq_len(burnin), ]
    }

    ## get vector of deviances
    ## using log likelihoods
    dev <- (-2) * sapply(seq_len(nrow(x)), function(y) f(as.numeric(x[y,argnames(f)])))

    ## estimate effective number of parameters
    dbar <- mean(dev)

    ## deviance of posterior means:
    ## evaluate deviance at the mean posterior estimate
    post.means <- sapply(argnames(f), function(y) mean(x[,y]))

    dhat <- -2 * f(post.means)

    pd <- dbar - dhat

    ## calculate dic
    dic <- dbar + pd
    unname(dic)
}


## function for making prior for bm
## just a wrapper for make.prior.uniform
make.prior.bm <- function(s2.lower, s2.upper){
    f <- make.prior.uniform(lower=s2.lower, upper=s2.upper)
    f
}



## define function to create a lognormal prior
## to be used in make.prior.ou
make.prior.lognormal <- function(ln.mean, ln.sd){
    function(x) sum(dlnorm(x, meanlog=ln.mean, sdlog = ln.sd, log=TRUE))
}



## function for making ou prior
## s2 is given a uniform prior
## alpha is given a lognormal prior
make.prior.ou <- function(s2.lower, s2.upper, ln.mean, ln.sd){
    p.s2 <- make.prior.uniform(lower=s2.lower, upper=s2.upper)
    p.al <- make.prior.lognormal(ln.mean=ln.mean, ln.sd=ln.sd)

    function(pars){
        p.s2(pars[1]) + p.al(pars[2])
    }
}



## function for making eb prior
## s2 is given a uniform prior
## a is givne a uniform prior
## the prior on s2 should be different than that of a
## s2 > 0 and a < 0
make.prior.eb <- function(s2.lower, s2.upper, a.lower, a.upper){
    p.s2 <- make.prior.uniform(lower=s2.lower, upper=s2.upper)
    p.a  <- make.prior.uniform(lower=a.lower, upper=a.upper)

    function(pars){
        p.s2(pars[1]) + p.a(pars[2])
    }
}


## create function for setting OU bounds for ML analyses
## this is the only one where bounds to be necessary
bnds.ou <- function(phy, states){
    ## index the terminal branches
    tips <- phy$edge[,2] <= Ntip(phy)
    ## find the shortest branch
    sht <- min(phy$edge.length[tips])
    ## upper bound for alpha
    ## half-life should be greater than shortest branch
    bnd.a <- log(2)/sht
    ## setting upper bound for sigsq at 2
    bnd <- c(2, bnd.a)
    bnd
}



## simple function for reading in angiosperm dataset
get.angio.data <- function(){
  lapply(dir("output/angio-data", full.names=TRUE), readRDS)
}


## data processing fxns
## function for reading in individual results
build.adequacy.results <- function(x){
    f <- read.csv(paste("output/temp", x, sep="/"), as.is=TRUE, row.names=1)
    d <- as.data.frame(t(f))
    rownames(d) <- NULL
    d
}


## simple function for unfactoring
## why is this not in base??
unfactor <- function(x)
    as.numeric(levels(x))[x]
