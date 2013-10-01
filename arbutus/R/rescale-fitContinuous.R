## fxns for rescaling tree according to fitContinuous type models
## adopted from rescale.phylo in GEIGER (written by Jon Eastman)


model.phylo.bm <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 | pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## rescale branch lengths according to sigsq
    sigsq <- pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= Ntip(phy)
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return tree
    phy
}





model.phylo.ou <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$alpha < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edge.height(phy)

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
    sigsq <- pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= Ntip(phy)
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}




## Rescale phylogeny according to early burst model
model.phylo.eb <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edge.height(phy)

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
model.phylo.lambda <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edge.height(phy)

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
    sigsq <- pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}





## kappa transformation
model.phylo.kappa <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$kappa < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## rescale branch lengths according to kappa
    kappa <- pars$kappa
    phy$edge.length <- phy$edge.length^kappa

    ## rescale branch lengths according to sigsq
    sigsq <- pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## add on SE to tip branches
    tips <- phy$edge[,2] <= Ntip(phy)
    phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE

    ## return phy
    phy
}





## delta transformation
model.phylo.delta <- function(phy, pars){
    ## check pars to make sure they are non-negative
    if (pars$sigsq < 0 || pars$delta < 0 || pars$SE < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    ht <- edge.height(phy)

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
    sigsq <- pars$sigsq

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

model.phylo.trend <- function(phy, pars)
    stop("trend model is not currently implemented")





## white noise transformation
## turn tree into star phylogeny
model.phylo.white <- function(phy, pars){
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







#' @method make.model.phylo fitC
#' @S3method make.model.phylo fitC
make.model.phylo.fitC <- function(x, ...){
    ## get model
    model <- x$type

    ## get tree
    phy <- x$data$phy

    ## get parameters
    pars <- x$pars

    ## use switch function to get appropriate transformation
    rphy <- switch(model,
                   BM = model.phylo.bm(phy, pars),
                   OU = model.phylo.ou(phy, pars),
                   EB = model.phylo.eb(phy, pars),
                   lambda = model.phylo.lambda(phy, pars),
                   kappa = model.phylo.kappa(phy, pars),
                   delta = model.phylo.delta(phy, pars),
                   trend = model.phylo.trend(phy, pars),
                   white = model.phylo.white(phy, pars),
               )

    ## return rescaled phylogeny
    rphy
}




#' @title Rescale phylogeny based on fitted model parameters
#'
#' @description Rescales phylogeny to form a 'unit.tree' using parameters
#' from fitted model
#'
#' @param x an object inherited from \code{\link{model.info}}
#'
#' @details This is a generic function which rescales the phylogeny based on the model
#' specific information. While the class and information may differ, the object must include
#' the 'phylo' object to be rescaled and model specific information such as model type and
#' parameter values.
#'
#' To include additional types of models, researchers will need to build a new \code{\link{model.info}}
#' function for the given model type. The output should be assigned a novel class.
#' A \code{make.model.phylo.classX} object will need to be created to perform the rescaling.
#'
#' @return a 'phylo' object that can be used to from a 'unit.tree' object
#'
#' @seealso \code{\link{as.unit.tree}}
#'
#' @export make.model.phylo
#'
#' @examples
#' ## finch data
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' data <- td$data[,"wingL"]
#'
#'
#' ## using just the given phylogeny
#' unit.tree.phy <- as.unit.tree(phy, data)
#'
#' ## fit Brownian motion model
#' ## using geiger's fitContinuous function
#' fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
#'                                  control=list(niter=10))
#'
#' ## get model info using internal arbutus function
#' info.bm <- arbutus:::model.info(fit.bm)
#'
#' ## rescale phylogeny based on model info
#' make.model.phylo(info.bm)
#' 
make.model.phylo <- function(x, ...)
    UseMethod("make.model.phylo")
