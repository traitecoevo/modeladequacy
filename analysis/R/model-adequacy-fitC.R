## Analysis script for looking at model adequacy across angiosperms
## Take 1 -- Oct.7 2013
##
## Note: eventually I want to turn this into a proper markdown file and
## include as a supplemental material

require(geiger)
require(arbutus)
source("extract_subtree_functions.R")


## Get bounds for OU
## For fitContinuous
bounds.ou <- function(phy, states){
    tips <- phy$edge[,2] <= Ntip(phy)
    sh.tip <- min(phy$edge.length[tips])

    ## quick and dirty estimate of sigsq
    unit.tree <- as.unit.tree(phy, data=states)
    sigsq <- sigsq.reml(unit.tree)

    ## upper bound of alpha parameter
    alpha.up <- 2 * sigsq / sh.tip

    bnd <- list(alpha=c(mn=0.0000001, mx=alpha.up))
    bnd
}


## Get bounds for EB
## for fitContinuous
bounds.eb <- function(phy){
    ht <- arbutus:::edge.height(phy)
    Tmax <- ht$start[Ntip(phy) + 1]

    ## lower bound for a parameter
    ## from Slater and Pennell 2013
    a.low <- log(10^(-5)) / Tmax

    bnd <- list(a=c(mn=a.low, mx=-0.0000001))
    bnd
}



## wrapper function which fits 3 fitContinuous models
## calculates aic weights
## applies model adequacy approach to best fit model
## returns aic weights and p values for summary stats
model.ad.angio.ml <- function(phy, states, SE){
    ## fit 3 models using fitContinuous
    fit.bm <- fitContinuous(phy, states, SE=SE, model="BM")
    fit.ou <- suppressWarnings(fitContinuous(phy, states,
                                             SE=SE, model="OU",
                                             bounds=bounds.ou(phy, states)))
    fit.eb <- suppressWarnings(fitContinuous(phy, states,
                                             SE=SE, model="EB",
                                             bounds=bounds.eb(phy)))

    ## get aic weights
    aic.bm <- fit.bm$opt$aic
    aic.ou <- fit.ou$opt$aic
    aic.eb <- fit.eb$opt$aic
    aic <- c(aic.bm, aic.ou, aic.eb)
    names(aic) <- c("BM", "OU", "EB")

    aic.w <- aicw(aic)

    ## pick the best model based on aic
    ## if there is a tie just pick one.
    best.fit <- sample(rownames(aic.w[which(aic.w[,"delta"] == 0),]), size=1)

    ## assess model adequacy
    pp <- switch(best.fit,
                 BM=phy.model.check(fit.bm),
                 OU=phy.model.check(fit.ou),
                 EB=phy.model.check(fit.eb))

    pval <- pval.summ.stats(pp)

    ## create output
    out <- c(aic.w["BM", "w"], aic.w["OU", "w"], aic.w["EB", "w"],
             best.fit, pval)
    tmp <- c("aic.w.bm", "aic.w.ou", "aic.w.eb", "model.used")
    names(out) <- c(tmp, names(pval))

    out
}
                      
 
