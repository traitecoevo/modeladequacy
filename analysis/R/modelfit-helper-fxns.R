library(arbutus)
library(diversitree)
library(multicore)

## function to calculate aic (or other weights) given a vector of aic scores
ic.weights <- function(x){
    
    d.x <- x - min(x)
    tmp <- exp(-0.5 * d.x)
    w   <- tmp / sum(tmp)
    names(w) <- names(x)
    w
}


## Define function for calculating dic from mcmcsamples
dic.mcmcsamples <- function(x, burnin=NULL){
    if (!inherits(x, "mcmcsamples"))
        stop("this function is only designed for diversitrees mcmcsamples object")

    f <- attr(x, "func")

    if (!is.null(burnin)){
        if (burnin == 0)
            stop("if burnin supplied, it cannot be 0")
        
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
    
## function for making prior for ou
## define function make.prior.lognormal
make.prior.lognormal <- function(ln.mean, ln.sd){
    function(x) sum(dlnorm(x, meanlog=ln.mean, sdlog = ln.sd, log=TRUE))
}


make.prior.ou <- function(s2.lower, s2.upper, ln.mean, ln.sd){
    p.s2 <- make.prior.uniform(lower=s2.lower, upper=s2.upper)
    p.al <- make.prior.lognormal(ln.mean=ln.mean, ln.sd=ln.sd)

    function(pars){
        p.s2(pars[1]) + p.al(pars[2])
    }
}


## make prior for early burst
make.prior.eb <- function(s2.lower, s2.upper, a.lower, a.upper){
    p.s2 <- make.prior.uniform(lower=s2.lower, upper=s2.upper)
    p.a  <- make.prior.uniform(lower=a.lower, upper=a.upper)

    function(pars){
        p.s2(pars[1]) + p.a(pars[2])
    }
}

    

## function for reading in data from source
get.angio.data <- function(){
    x <- readRDS(file="output/angio-trait-data-all.rds")
    x
}

    
    
