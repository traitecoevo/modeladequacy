library(arbutus)
library(diversitree)
library(nlme)

## simulate some trait univariate data under OU along a phylogeny

tree <- sim.bdtree(stop="taxa", n=100)

ou.data <- sim.char(rescale(tree, model="OU", alpha=1), par=1, model="BM")[,,]


## fit with BM

## using fitContinuous
fit.C <- fitContinuous(phy=tree, dat=ou.data, SE=0, model="BM")

## using diversitree
bmlik <- make.bm(tree=tree, states=ou.data)

fit.D <- find.mle(bmlik, x.init=1)


## now perform a model adequacy test using default parameters
## under the hood:
## 1. parse model fit from fitted object
## 2. rescale tree according to model and parameters
## 3. compute summ stats on observed data
## 4. simulate 1000 datasets
## 5. calc summary stats on simulated data
## 6. compare simulated to observed summary stats

## take a look at the main fxn phyModelCheck
phyModelCheck

## you will notice that there are two fxns modelinfo() and modelphylo()
## which are generic. To add more additional models all one needs to do is write a new
## modelinfo.whatever() to pull out the info from a object of class whatever
## and a modelphylo.whatever() fxn to tell it how to rescale the phylogeny based on the parameters


## Currently only outputting p-vlaues but I could change this


## using geiger fitted object
res.fC <- phyModelCheck(fit.C)

res.fC

## using diversitree fitted object
## due to diversitree's internals, you have to also supply likelihood function
res.dt <- phyModelCheck(fit.D, lik=bmlik)

res.dt


## or if you want to just rescale the tree yourself and not use fitted object
my.tree <- rescale(tree, model="BM", sigsq=fit.C$opt$sigsq)

## and include data in phyModelCheck
res.mytree <- phyModelCheck(my.tree, data=ou.data)

res.mytree



## amd say, you were only interested in calculating the mean value of the pics
## you can just create a fxn
my.fxn <- function(x){
    mean(x$pics[,"contrasts"])
}

## and then include it as a list (currently need to do this the long way, but will fix it)
ut <- as.unit.tree(fit.C)
ss.obs <- summStat(ut, stats=list(mean=my.fxn))
sims <- sim.charUnit(ut, nsim=1000)
ss.sim <- summStat(sims, stats=list(mean=my.fxn))
p.value <- compare.summStat(ss.obs, ss.sim)








## Now is where it gets really good
## Simulate two characters
set.seed(1)
x <- sim.char(tree, par=1, model="BM")[,,]
y <- sim.char(tree, par=1, model="BM")[,,]

dat <- cbind.data.frame(x,y)


## now we can do PGLS analysis

## using Brownian Motion
pgls.bm <- gls(x~y, data=dat, correlation=corBrownian(phy=tree))

## using a Lambda transform
pgls.lamb <- gls(x~y, data=dat, correlation=corPagel(value=1, phy=tree))


## and use same fxns as before
## pulls out the tree, residuals and parameters of correlation structure

res.pglsbm <- phyModelCheck(pgls.bm)

res.pglsbm


## and again with the lambda transform

res.pglslamb <- phyModelCheck(pgls.lamb)

res.pglslamb
