## Analysis script for looking at model adequacy across angiosperms
## Take 1 -- Oct.7 2013
##
## Note: eventually I want to turn this into a proper markdown file and
## include as a supplemental material

require(geiger)
require(arbutus)
require(multicore)
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

    ## create output
    out <- c(aic.w["BM", "w"], aic.w["OU", "w"], aic.w["EB", "w"],
             best.fit, pval)
    tmp <- c("aic.w.bm", "aic.w.ou", "aic.w.eb", "model.used")
    names(out) <- c(tmp, names(pval))

    out
}







## READ IN THE DATA


## read in bigtree
t <- read.tree(file="tempo_scrubbed_CONSTRAINT_rooted.dated.tre")

## ANALYSIS OF SLA DATA
sla.raw <- read.csv(file="species_mean_sla.csv")
sla <- sla.raw[,"x"]

## log the sla using base 10
## same units SE was calculated in
sla <- log10(sla)
names(sla) <- sla.raw[,"X"]

## SE
sla.se <- 0.1024202








### TIME SLICE ANALYSIS


slice.and.model.select<-function(age,smaller.tree=smaller.tree,sr.min,trait.vec=l.sla){
  tree.list<-time.slice.tree(age,smaller.tree,sr.min)
  ms.out<-lapply(tree.list,modelad.ml, states=trait.vec, SE=0.1024202)
  return(ms.out)
} 


l.sla<-log10(read.csv("../output/species_mean_sla.csv",row.names=1))
tree<-read.tree("../data/tempo_scrubbed_CONSTRAINT_rooted.dated.tre")
smaller.tree<-geiger:::.drop.tip(tree,tree$tip.label[!tree$tip.label%in%row.names(l.sla)])

#age of full tree is 400.7877 million years, and time slices measures from the root toward the tips
#time.slices<-c(0.7877,50.7877,100.7877,150.7877,200.7877,250.7877,300.7877,350.7877,375.7877)
time.slices<-c(0.7877,250.7877,350.7877)
diff.out<-lapply(X=time.slices,FUN=slice.and.model.select,smaller.tree=smaller.tree,sr.min=50)








## function which pulls out data sets by clade from
## full tree and dataset
## assume rownames of dataset are taxon labels
## rank can be 'family' or 'order' (will add genus later)
## min.size is the minimum clade size
treedata.taxon <- function(phy, data, rank="family", min.size=20){
    ## drop tips that are not in data set
    ## to avoid dropping same tips over and over
    tt <- phy$tip.label[-which(phy$tip.label %in% names(data))]
    phy <- geiger:::.drop.tip(phy, tip=tt)

    if (rank == "family"){
        ## get all family nodes in tree
        tax <- phy$node.label[grep("[A-z]+ceae", phy$node.label, perl=TRUE)]
        ## extract subtrees
        trees <- lapply(tax, function(x) extract.clade(phy, node=x))
        ## use treedata to match to family level
        td <- lapply(trees, function(x) treedata(phy=x, data=data))
        names(td) <- tax
        ## get number of tips
        tips <- lapply(td, function(x) Ntip(x$phy))
        ## extract those which meet threshold
        dd <- td[tips >= min.size]
    }
    if (rank == "order"){
        ## get all ordinal nodes in tree
        tax <- phy$node.label[grep("[A-z]+ales", phy$node.label, perl=TRUE)]
        ## extract subtrees
        trees <- lapply(tax, function(x) extract.clade(phy, node=x))
        ## use treedata to match to family level
        td <- lapply(trees, function(x) treedata(phy=x, data=data))
        names(td) <- tax
        ## get number of tips
        tips <- lapply(td, function(x) Ntip(x$phy))
        ## extract those which meet threshold
        dd <- td[tips >= min.size]
    }

    dd
}












## get family data
fam.data <- treedata.taxon(phy=t, data=sla, rank="family", min.size = 20)

## get order data
ord.data <- treedata.taxon(phy=t, data=sla, rank="order", min.size = 20)


## family level model adequacy

res.fam <- lapply(fam.data, function(x) modelad.ml(x$phy, x$data, SE = sla.se))
res.fam <- do.call(rbind, res.fam)
trait <- rep("sla", nrow(res.fam))
res.fam <- cbind.data.frame(trait, res.fam)
write.csv(res.fam, file="sla-family-ml-res.csv")

res.ord <- lapply(ord.data, function(x) modelad.ml(x$phy, x$data, SE = sla.se))
res.ord <- do.call(rbind, res.ord)
trait <- rep("sla", nrow(res.ord))
res.ord <- cbind.data.frame(trait, res.ord)
write.csv(res.ord, file="sla-order-ml-res.csv")







