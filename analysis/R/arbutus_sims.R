## load in model adequacy fitting fxns
source("R/model-adequacy-fit.R")
source("R/sim-helper.R")

## create directory: note for reproducibility, this should be cleaned up
dir.create(path.sim())

mc.cores <- 24

## Sims start here
n.taxa <- c(50, 100, 200)

## SE
se <- c(0, 0.05)

## Just a small sampling of parameter space
## BM sims
pars.bm.1 <- list(model="BM", sigsq=1, alpha=0, a=0)
pars.ou.1 <- list(model="OU", sigsq=1, alpha=1, a=0)
pars.ou.2 <- list(model="OU", sigsq=1, alpha=2, a=0)
pars.ou.3 <- list(model="OU", sigsq=1, alpha=4, a=0)
pars.eb.1 <- list(model="EB", sigsq=1, alpha=0, a=log(0.01))
pars.eb.2 <- list(model="EB", sigsq=1, alpha=0, a=log(0.02))
pars.eb.3 <- list(model="EB", sigsq=1, alpha=0, a=log(0.04))

all.pars <- list(pars.bm.1, pars.ou.1, pars.ou.2, pars.ou.3,
                 pars.eb.1, pars.eb.2, pars.eb.3)

for (j in seq_len(length(all.pars))){
    for (i in seq_len(length(se))){
        for (k in seq_len(length(n.taxa))){
            pars <- all.pars[[j]]
            pars$SE <- se[i]
            model <- pars$model
            filename <- paste(model, pars$sigsq, pars$alpha, pars$a, pars$SE)
            model_ad_sim_multi(pars, model, "ml", n.taxa[k], 1000, write=TRUE, filename)
        }
    }
}

