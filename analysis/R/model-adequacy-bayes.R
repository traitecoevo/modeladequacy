## Script for performing model adequacy using a mcmc
library(arbutus)
library(LaplacesDemon)
#source("R/extract_subtree_functions.R")
#source("R/read-data-functions.R")
library(diversitree)

## function for calculating Deviance Information Criterion
## from mcmcsamples

dic.mcmcsamples <- function(x, lik){
    ## get vector of deviances
    dev <- -2 * x$p

    ## estimate effective number of parameters
    ## posterior mean of the deviance
    dbar <- mean(dev)

    ## deviance of posterior means
    ## evaluate deviance at the mean posterior estimates
    post.means <- sapply(argnames(lik), function(y) mean(x[,y]))
    
    dhat <- -2 * lik(post.means)

    pd <- dbar - dhat

    ## calculate dic
    dic <- dbar + pd

    unname(dic)
}









## function for making a half-cauchy distribution
## using LaplacesDemon pkg
make.prior.halfcauchy <- function(scale, log=TRUE){
    function(x) sum(dhalfcauchy(x, scale=scale, log=log))
}



## function for making ou priors
## mixture of uniform and half-cauchy
make.prior.ou <- function(s2.low, s2.high, scale, theta.low, theta.high){
    pr.s2 <- make.prior.uniform(s2.low, s2.high)
    pr.alpha <- make.prior.halfcauchy(scale)
    pr.theta <- make.prior.uniform(theta.low, theta.high)

    function(pars){
        pr.s2(pars[1]) + pr.alpha(pars[2]) + pr.theta(pars[3])
    }
}




## function for making eb priors
make.prior.eb <- function(s2.low, s2.high, a.low, a.high){
    pr.s2 <- make.prior.uniform(s2.low, s2.high)
    pr.a <- make.prior.uniform(a.low, a.high)

    function(pars){
        pr.s2(pars[1]) + pr.a(pars[2])
    }
}



## eventually to be wrapped for data analysis

## for now just simulate data
phy <- tree.bd(pars=c(1,0), max.taxa=200)
states <- sim.character(t, 1, x0=0)


## BM
## create likelihood fxn

lik.bm <- make.bm(phy, states, control=list(method="pruning"))

## bm prior
prior.bm <- make.prior.uniform(lower=0, upper=10)

## use find.mle to get reasonable starting values
fit.bm <- find.mle(lik.bm, x.init = 1)

## run temporary chain to get a better w
tmp.bm <- mcmc(lik.bm, x.init = fit.bm$par, nsteps = 100, prior=prior.bm,
               lower=0, w=1, print.every=0)
w.bm <- diff(sapply(tmp.bm[2], range))

samples.bm <- mcmc(lik.bm, fit.bm$par, w=w.bm, nsteps = 10000, prior=prior.bm,
                   print.every=0)




## OU
## create likelihood fxn

lik.ou <- make.ou(phy, states, control=list(method="pruning"))

## ou priors
## half cauchy with scale of 25 is used as a standard,
## weakly informative scale parameter prior
## can fiddle with this
prior.ou <- make.prior.ou(s2.low=0, s2.high=10, scale=25,
                          theta.low = min(states), theta.high = max(states))

## use find.mle to get reasonable starting values
fit.ou <- find.mle(lik.ou, x.init = c(1,0.1,mean(states)))

## run temporary chain to get a better w
tmp.ou <- mcmc(lik.ou, fit.ou$par, nsteps = 100, prior=prior.ou,
               w=1, print.every=0)

w.ou <- diff(sapply(tmp.ou[2:4], range))

samples.ou <- mcmc(lik.ou, fit.ou$par, nsteps = 10000, prior=prior.ou,
                   w=w.ou, print.every=0)



## EB
## create likelihood fxn

lik.eb <- make.eb(phy, states, control=list(method="pruning"))

## eb priors.ou
tmax <- max(branching.times(phy))
prior.eb <- make.prior.eb(s2.low=0, s2.high=10, a.low=log(10^(-5))/tmax, a.high=0)

## use find.mle to get reasonable starting values
fit.eb <- find.mle(lik.eb, x.init = c(1, -0.1))

## run temporary chain to get better w
tmp.eb <- mcmc(lik.eb, fit.eb$par, nsteps = 100, prior=prior.eb,
               w=1, print.every=0)

w.eb <- diff(sapply(tmp.eb[2:3], range))

samples.eb <- mcmc(lik.eb, fit.eb$par, nsteps = 10000, prior=prior.eb,
                   w=w.eb, print.every=0)





## remove burnin before making model comparisons
burn <- 1000
res.bm <- samples.bm[-seq_len(burn),]
res.ou <- samples.ou[-seq_len(burn),]
res.eb <- samples.eb[-seq_len(burn),]




## get DIC for all three
dic.bm <- dic.mcmcsamples(res.bm, lik.bm)
dic.ou <- dic.mcmcsamples(res.ou, lik.ou)
dic.eb <- dic.mcmcsamples(res.eb, lik.eb)
dic <- c(dic.bm, dic.ou, dic.eb)
names(dic) <- c("BM", "OU", "EB")

## get best model (min dic)
## if tie, sample
best.fit <- sample(names(dic[which(dic == min(dic))]), size=1)


## run model adequacy using best fit model
pp <- switch(best.fit,
             BM=phy.model.check(res.bm, lik=lik.bm, n.samples=1000,
                 check=FALSE),
             OU=phy.model.check(res.ou, lik=lik.ou, n.samples=1000,
                 check=FALSE),
             EB=phy.model.check(res.eb, lik=lik.eb, n.samples=1000,
                 check=FALSE))


## get pvalues
pval <- pval.summ.stats(pp)

## get mean mahalanobis distance
m.dist <- mean(mv.summ.stats(pp))

## create output
out <- c(dic.bm, dic.ou, dic.eb, best.fit, m.dist, pval)
nam <- c("dic.bm", "dic.ou", "dic.eb", "model.used", "mv.modelad")
names(out) <- c(nam, names(pval))
out 






