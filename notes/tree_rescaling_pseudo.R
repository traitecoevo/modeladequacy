## This is psuedocode for the general structure of the package
## Aug 29, 2013

## Matt Pennell

## The goal here is to make this user friendly, flexible and robust
## at the same time

## here is what i am thinking


## 1. user fits model using whatever software they like
## will use fitContinuous here

fit <- fitContinuous(phy=phy, dat=dat, model="BM")


## 2. create parsing functions which take output of model fit from program X

res <- parse.fitContinuous(fit)

## the parser functions reads the output of the model fitting and returns a list of rescaling functions
## if single model, such as BM, returns a single fxn. if more complex model, returns a list equal to the number of branches

## so if the model was originally BM:
## res would be return
[[1]]
edge.BM(sigsq)

## each function returns will inturn be used to create fxns which take a branch length (as well as total tree depth for OU, lambda, etc. models)
foo <- edge.BM(sigsq)

## foo() returns
fxn(branch.start, branch.end)




## The edge fxns can then be used to rescale the tree when building hte unit.tree

## The as.unit.tree fxn does the following:
## If fxns=NULL

unit.tree <- as.unit.tree(phy=phy, data=data, fxns=NULL){
    if (!inherits(x, "phylo"))
        stop("phylogney must be of class phylo")

    if (fxns=NULL){ ## use the tree as given. Do not rescale
        
        ## check tree and data to make sure they match
        ## calculate pics
        ## create unit.tree object consiting of $phy, $data, $pics

    } else {

        ## get start times and end times for every branch
        ## get total tree depth

        ## use internal function which takes the fxns and creates a scaling fxn for everybranch in the tree
        ## if only one fxn supplied, apply ot all branches
        
        ## scale all branches according to specified branch wise fxns

        ## calculate pics on new tree
        ## create unit.tree object consitsting of $phy, $data, $pics
        
        
    }

}


unit.tree <- as.unit.tree(phy, data, fxns=edgeBM(sigsq))


## calculate summary stats on using your observed data

ss.obs <- summStat(unit.tree, stats=NULL) ## use default stats. User supplied stats can also be used

## simulate n data sets on unit tree
sims <- sim.charUnit(unit.tree, n=1000)

## calculate summary stats on simulated data sets

ss.sim <- summStat(sims, stats=NULL)

## compare the observed summary statistics with the simulated summary statistics

p.value <- compare.summStat(ss.obs, ss.sim)


