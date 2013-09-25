## fxns for rescaling tree according to fitContinuous type models
## takes modelinfo objects

## model specific rescaling fxns
## adopted from rescale.phylo in GEIGER (written by Jon Eastman)

## import heights.phylo and .paths.phylo from GEIGER








modelphyloBM <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 | pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## rescale branch lengths according to sigsq
    sigsq = pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= Ntip(phy)
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return tree
    phy
}





modelphyloOU <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$alpha < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edgeHeight(phy)

    ## get tree depth
    ## same as max(branching.times(phy)) but works for non-ultrametric trees
    N <- Ntip(phy)
    Tmax <- ht$start[N+1]

    ## get braching times in terms of tree depth
    idx <- match(1:nrow(ht), phy$edge[,2])
    t1 <- Tmax - ht$end[phy$edge[idx,1]]
    t2 <- ht$start - ht$end + t1

    ## rescale branch lengths according to alpha
    alpha <- pars$alpha

    bl <- (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t2)) * (1 - exp(-2 * alpha * t2)) -
          (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t1)) * (1 - exp(-2 * alpha * t1))

    phy$edge.length <- bl[phy$edge[,2]]

    ## rescale branch lengths according to sigsq
    sigsq = pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= Ntip(phy)
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}




## Rescale phylogeny according to early burst model
modelphyloEB <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edgeHeight(phy)

    ## get tree depth
    ## same as max(branching.times(phy)) but works for non-ultrametric trees
    N <- Ntip(phy)
    Tmax <- ht$start[N+1]

    ## get braching times in terms of tree depth
    idx <- match(1:nrow(ht), phy$edge[,2])
    t1 <- Tmax - ht$end[phy$edge[idx,1]]
    t2 <- ht$start - ht$end + t1

    ## rescale tree according to a parameter
    a <- pars$a

    if (a != 0){
        bl <- (exp(a * t2) - exp(a * t1)) / a
        phy$edge.length <- bl[phy$edge[,2]]
    }

    ## rescale branch lengths according to sigsq
    sigsq <- pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= Ntip(phy)
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}





## lambda transformation
modelphyloLambda <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edgeHeight(phy)

    ## get tree depth
    ## same as max(branching.times(phy)) but works for non-ultrametric trees
    N <- Ntip(phy)
    Tmax <- ht$start[N+1]

    ## index edge matrix by tips
    tips <- match(1:N, phy$edge[,2])

    ## get all path lengths
    ## call internal geiger fxn
    pp <- geiger:::.paths.phylo(phy)

    ## rescale all branches according to lambda value
    lambda <- pars$lambda

    bl <- phy$edge.length*lambda

    ## readjust the tip branches
    bl[tips] <- bl[tips] + (pp - (pp * lambda))
    phy$edge.length <- bl

    ## check to make sure no negative branch lengths produced
    if (any(phy$edge.length < 0))
        warning("negative branch lengths in tranformed tree. Lambda may be too large")

    ## rescale branch lengths according to sigsq
    sigsq = pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}





## kappa transformation
modelphyloKappa <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$kappa < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## rescale branch lengths according to kappa
    kappa <- pars$kappa
    phy$edge.length <- phy$edge.length^kappa

    ## rescale branch lengths according to sigsq
    sigsq = pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= Ntip(phy)
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}





## delta transformation
modelphyloDelta <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$delta < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edgeHeight(phy)

    ## get tree depth
    ## same as max(branching.times(phy)) but works for non-ultrametric trees
    N <- Ntip(phy)
    Tmax <- ht$start[N+1]

    ## get braching times in terms of tree depth
    ## note: Geiger's delta transformation does not index the branches (like OU)
    ## need to check if this is correct
    t <- Tmax - ht$start
    l <- ht$start - ht$end

    ## rescale the branches according to delta
    delta <- pars$delta

    bl <- (t + l)^delta - t^delta
    phy$edge.length <- bl[phy$edge[,2]]

    ## rescale branch lengths according to sigsq
    sigsq = pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= N
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}




## trend transformation
## NEED TO WRITE THIS FXN LATER
## placeholder fxn

modelphyloTrend <- function(phy, pars)
    stop("trend model is not currently implemented")





## white noise transformation
## turn tree into star phylogeny
modelphyloWhite <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0)
        stop("Parameters need to be non-negative")

    ## collapse to star phylogeny
    tips <- Ntip(phy)
    phy$edge.length[] <- 0
    phy$edge.length[phy$edge[,2] <= tips] <- 1

    ## rescale according to sigsq
    sigsq <- pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## return phy
    phy
}





## modelPhylo fxn for rescaling based on fitContinuous styles model.fitC
## includes fitContinuous, diversitree, gls, pgls, MCMCglmm, fitContinuousMCMC, etc.

modelphylo.model.fitC <- function(x, ...){
    ## get model
    model <- x$type

    ## get tree
    phy <- x$data$phy

    ## get parameters
    pars <- x$pars

    ## use switch function to get appropriate transformation
    rphy <- switch(model,
                   BM = modelphyloBM(phy, pars),
                   OU = modelphyloOU(phy, pars),
                   EB = modelphyloEB(phy, pars),
                   lambda = modelphyloLambda(phy, pars),
                   kappa = modelphyloKappa(phy, pars),
                   delta = modelphyloDelta(phy, pars),
                   trend = modelphyloTrend(phy, pars),
                   white = modelphyloWhite(phy),
               )

    ## return rescaled phylogeny
    rphy
}
