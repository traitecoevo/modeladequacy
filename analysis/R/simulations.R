## load in model adequacy fitting fxns
source("R/model-adequacy-fit.R")
source("R/sim-helper.R")


## create directory: note for reproducibility, this should be cleaned up
dir.create(path.sim())

n.cores <- 20

## Sims start here
n.taxa <- c(50, 100, 200)

## Just a small sampling of parameter space
## Not intended to be comprehensive

## BM sims
pars.bm.nose <- list(sigsq=1, alpha=0, a=0, SE=0)
pars.bm.se <- list(sigsq=1, alpha=0, a=0, SE=0.05)
lapply(n.taxa, function(x) sim_ad(pars.bm.nose, "BM", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.bm.se, "BM", "ml", x, 1000))


## OU sims
## 3 parameter sets
pars.ou.nose.1 <- list(sigsq=1, alpha=1, a=0, SE=0)
pars.ou.nose.2 <- list(sigsq=1, alpha=2, a=0, SE=0)
pars.ou.nose.3 <- list(siqsq=1, alpha=4, a=0, SE=0)
pars.ou.se.1 <- list(sigsq=1, alpha=1, a=0, SE=0.05)
pars.ou.se.2 <- list(sigsq=1, alpha=2, a=0, SE=0.05)
pars.ou.se.3 <- list(siqsq=1, alpha=4, a=0, SE=0.05)
lapply(n.taxa, function(x) sim_ad(pars.ou.nose.1, "OU", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.ou.nose.2, "OU", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.ou.nose.3, "OU", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.ou.se.1, "OU", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.ou.se.2, "OU", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.ou.se.3, "OU", "ml", x, 1000))

## EB sims
## 3 parameter sets
pars.eb.nose.1 <- list(sigsq=1, alpha=0, a=log(0.01), SE=0)
pars.eb.nose.1 <- list(sigsq=1, alpha=0, a=log(0.02), SE=0)
pars.eb.nose.1 <- list(sigsq=1, alpha=0, a=log(0.04), SE=0)
pars.eb.nose.1 <- list(sigsq=1, alpha=0, a=log(0.01), SE=0.05)
pars.eb.nose.1 <- list(sigsq=1, alpha=0, a=log(0.02), SE=0.05)
pars.eb.nose.1 <- list(sigsq=1, alpha=0, a=log(0.04), SE=0.05)
lapply(n.taxa, function(x) sim_ad(pars.eb.nose.1, "EB", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.eb.nose.2, "EB", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.eb.nose.3, "EB", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.eb.se.1, "EB", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.eb.se.2, "EB", "ml", x, 1000))
lapply(n.taxa, function(x) sim_ad(pars.eb.se.3, "EB", "ml", x, 1000))







        






