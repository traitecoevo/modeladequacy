## Script for performing model adequacy using a mcmc
library(arbutus)
library(LaplacesDemon)
source("R/extract_subtree_functions.R")
source("R/read-data-functions.R")
library(diversitree)
library(multicore)

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
#phy <- tree.bd(pars=c(1,0), max.taxa=200)
# states <- sim.character(t, 1, x0=0)



modelad.bayes <- function(phy, states, SE){

## BM
## create likelihood fxn

lik.bm <- make.bm(phy, states, states.sd = SE, control=list(method="pruning", backend="C"))

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

lik.ou <- make.ou(phy, states, states.sd = SE, control=list(method="pruning", backend="C"))

## ou priors
## half cauchy with scale of 25 is used as a standard,
## weakly informative scale parameter prior
## can fiddle with this
prior.ou <- make.prior.ou(s2.low=0, s2.high=10, scale=25,
                          theta.low =0, theta.high = max(states))

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

lik.eb <- make.eb(phy, states, states.sd = SE, control=list(method="pruning", backend="C"))

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

}




modelad.bayes.clade <- function(td, se, rank, trait){
    states <- td$data[,1]
    ## perform model adequacy using bayesian analysis
    tmp <-  modelad.bayes(phy=td$phy, states=states, SE=se)

    ## get age of taxa
    age <- arbutus:::edge.height(td$phy)$start[Ntip(td$phy)+1]

    size <- Ntip(td$phy)

    out <- c(rank, trait, size, age, tmp)
    names(out) <- c("rank", "trait", "size", "age", names(tmp))

    out
}



# modelad.bayes.slice <- function(td, trait.name){
    ## perform model adequacy using bayesian analysis
  #  tmp <- modelad.bayes(phy=x$phy, states=x$data, SE=tree.states$SE)

    ## get age of taxa
  #  age <- sapply(td, function(x) return(arbutus:::edge.height(x$phy)$start[Ntip(x$phy)+1]))

  #  taxa <- rank <- rep(NA, length(age))

  #  trait <- rep(trait.name, length(age))

  #  size <- sapply(td, function(x) Ntip(x$phy))

  #  out <- cbind.data.frame(taxa, rank, trait, size, age, tmp)
  #  colnames(out) <- c("taxa", "rank", "trait", "size", "age", colnames(tmp))

  #  out
    
# }



get.timeslice.treedata <- function(phy, states, age, min.size){
     ## get list of sliced trees
    trees <- time.slice.tree(time.slice=age, temp.tree=phy, sr=min.size)

    ## append correct data to each tree
    ## do this so bounds are properly set within modelad.ml
    td <- lapply(trees, function(x) treedata(phy=x, data=states))

    td
}





## detach geiger and coda to avoid conflict with mcmc

# time.slices <- c(0.2697, 50.2697, 100.2697, 150.2697, 200.2697)

## SLA
## read in data
sla <- get.sla.data()



td.sla.fam <- treedata.taxon(phy=sla$phy, data=sla$states,
                             rank="family", min.size = 20)
sla.fam <- mclapply(td.sla.fam, function(x) modelad.bayes.clade(x, se=sla$SE, rank="family", trait="sla"),
                    mc.preschedule = FALSE, mc.cores = 3)

sla.fam <- do.call(rbind, sla.fam)
rownames(sla.fam) <- names(td.sla.fam)
sla.fam <- as.data.frame(sla.fam)

write.csv(sla.fam, file="output/results-bayes-angio-sla-family.csv")





td.sla.ord <- treedata.taxon(phy=sla$phy, data=sla$states,
                             rank="order", min.size = 20)
sla.ord <- mclapply(td.sla.ord, function(x) modelad.bayes.clade(x, se=sla$SE, rank="order", trait="sla"),
                    mc.preschedule = FALSE, mc.cores = 3)

sla.ord <- do.call(rbind, sla.ord)
rownames(sla.ord) <- names(td.sla.ord)
sla.ord <- as.data.frame(sla.ord)

write.csv(sla.ord, file="output/results-bayes-angio-sla-order.csv")







## seedmass

sm <- get.seedmass.data()



td.sm.fam <- treedata.taxon(phy=sm$phy, data=sm$states,
                             rank="family", min.size = 20)
sm.fam <- mclapply(td.sm.fam, function(x) modelad.bayes.clade(x, se=sm$SE, rank="family", trait="seedMass"),
                   mc.preschedule = FALSE, mc.cores = 3)

sm.fam <- do.call(rbind, sm.fam)
rownames(sm.fam) <- names(td.sm.fam)
sm.fam <- as.data.frame(sm.fam)

write.csv(sm.fam, file="output/results-bayes-angio-seedmass-family.csv")




td.sm.ord <- treedata.taxon(phy=sm$phy, data=sm$states,
                             rank="order", min.size = 20)
sm.ord <- mclapply(td.sm.ord, function(x) modelad.bayes.clade(x, se=sm$SE, rank="order", trait="seedMass"),
                   mc.preschedule = FALSE, mc.cores = 3)

sm.ord <- do.call(rbind, sm.ord)
rownames(sm.ord) <- names(td.sm.ord)
sm.ord <- as.data.frame(sm.ord)

write.csv(sm.ord, file="output/results-bayes-angio-seedmass-order.csv")






## leafN

ln <- get.leafn.data()



td.ln.fam <- treedata.taxon(phy=ln$phy, data=ln$states,
                             rank="family", min.size = 20)
ln.fam <- mclapply(td.ln.fam, function(x) modelad.bayes.clade(x, se=ln$SE, rank="family", trait="leafN"),
                   mc.preschedule = FALSE, mc.cores = 24)

ln.fam <- do.call(rbind, ln.fam)
rownames(ln.fam) <- names(td.ln.fam)
ln.fam <- as.data.frame(ln.fam)

write.csv(ln.fam, file="output/results-bayes-angio-leafn-family.csv")




td.ln.ord <- treedata.taxon(phy=ln$phy, data=ln$states,
                             rank="order", min.size = 20)
ln.ord <- mclapply(td.ln.ord, function(x) modelad.bayes.clade(x, se=ln$SE, rank="order", trait="leafN"),
                   mc.preschedule = FALSE, mc.cores = 24)

ln.ord <- do.call(rbind, ln.ord)
rownames(ln.ord) <- names(td.ln.ord)
ln.ord <- as.data.frame(ln.ord)

write.csv(ln.ord, file="output/results-bayes-angio-leafn-order.csv")










