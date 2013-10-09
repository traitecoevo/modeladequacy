## Analysis script for looking at model adequacy across angiosperms
## Take 1 -- Oct.7 2013
##
## Note: eventually I want to turn this into a proper markdown file and
## include as a supplemental material

require(geiger)
require(arbutus)
cd <- getwd()
esf <- file.path(cd, "R", "extract_subtree_functions.R")
rdf <- file.path(cd, "R", "read-data-functions.R")
source(esf)
source(rdf)


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




## function for logging ks values
log.ks <- function(x){
    x$summ.stats.obs[,"ks.dstat"] <- log(x$summ.stats.obs[,"ks.dstat"])
    x$summ.stats.sim[,"ks.dstat"] <- log(x$summ.stats.sim[,"ks.dstat"])
    x
}


## wrapper function which fits 3 fitContinuous models
## calculates aic weights
## applies model adequacy approach to best fit model
## returns aic weights and p values for summary stats
modelad.ml <- function(phy, states, SE){
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

    ss <- log.ks(pp)

    obs <- as.matrix(ss$summ.stats.obs)
    sim <- as.matrix(ss$summ.stats.sim)
    cv.sim <- cov(sim)

    m <- mahalanobis(x=obs, center = colMeans(sim), cov=cv.sim)

    ## create output
    out <- c(aic.w["BM", "w"], aic.w["OU", "w"], aic.w["EB", "w"],
             best.fit, m, pval)
    tmp <- c("aic.w.bm", "aic.w.ou", "aic.w.eb", "model.used", "mv.modelad")
    names(out) <- c(tmp, names(pval))

    out
}







## function for doing time slice analyses

modelad.ml.slice <- function(tree.states, age, sr.min, trait.name){
    ## get list of sliced trees
    trees <- time.slice.tree(time.slice=age, temp.tree=tree.states$phy, sr=sr.min)

    ## append correct data to each tree
    ## do this so bounds are properly set within modelad.ml
    td <- lapply(trees, function(x) treedata(phy=x, data=tree.states$states))

    ## perform model adequacy using ml
    res <- lapply(td, function(x) modelad.ml(phy=x$phy, states=x$data, SE=tree.states$SE))
    tmp <- do.call(rbind, res)

    ## get age of taxa
    age <- sapply(td, function(x) return(arbutus:::edge.height(x$phy)$start[Ntip(x$phy)+1]))

    taxa <- rank <- rep(NA, length(age))

    trait <- rep(trait.name, length(age))

    size <- sapply(td, function(x) Ntip(x$phy))

    out <- cbind.data.frame(taxa, rank, trait, size, age, tmp)
    colnames(out) <- c("taxa", "rank", "trait", "size", "age", colnames(tmp))

    out
    
}







## function for doing cladewise analyses

modelad.ml.clade <- function(tree.states, rank, min.size, trait.name){
    ## get list of clade trees and data
    td <- treedata.taxon(phy=tree.states$phy, data=tree.states$states,
                         rank=rank, min.size=min.size)

    ## perform model adequacy using ml
    res <- lapply(td, function(x) modelad.ml(phy=x$phy, states=x$data, SE=tree.states$SE))
    tmp <- do.call(rbind, res)

    ## get age of taxa
    age <- sapply(td, function(x) return(arbutus:::edge.height(x$phy)$start[Ntip(x$phy)+1]))

    taxa <- names(td)

    rank <- rep(rank, length(age))

    trait <- rep(trait.name, length(age))

    size <- sapply(td, function(x) Ntip(x$phy))

    out <- cbind.data.frame(taxa, rank, trait, size, age, tmp)
    colnames(out) <- c("taxa", "rank", "trait", "size", "age", colnames(tmp))
    rownames(out) <- NULL

    out
}

    

    


## ANALYSIS STARTS HERE

## SLA
## read in data
sla <- get.sla.data(file.path(cd, "data", "tempo_scrubbed_CONSTRAINT_rooted.dated.tre"),
                           file.path(cd, "data", "species_mean_sla.csv"))


## clade analyses
fam <- modelad.ml.clade(sla, rank="family", min.size=20, trait.name="sla")
write.csv(fam, file=file.path(cd, "output", "results-ml-sla-family.csv"))

ord <- modelad.ml.clade(sla, rank="order", min.size=20, trait.name="sla")
write.csv(ord, file=file.path(cd, "output", "results-ml-sla-order.csv"))


## time slice analysis
time.slices <- c(0.7877, 50.7877, 100.7877, 150.7877, 200.7877,
                 250.7877, 300.7877, 350.7877, 375.7877)
ts <- lapply(time.slices, function(x) modelad.ml.slice(sla, age=x,
                                                       sr.min=20, trait.name="sla"))
ts.res <- do.call(rbind, ts)
write.csv(ts.res, file=file.path(cd, "output", "results-ml-sla-timeslice.csv"))
