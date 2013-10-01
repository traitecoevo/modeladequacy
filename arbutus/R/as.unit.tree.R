#' @title Create unit.tree object
#'
#' @description Generic method for creating a 'unit.tree' object
#'
#' @param x fitted model object or 'phylo' object
#'
#' @param ... additional arguments
#'
#' @details This function is a generic function which takes a fitted
#' model object (from a number of different packages), rescales the phylogeny
#' based on the model fitted parameters, computes the contrasts on the rescaled
#' phylogeny and returns an object of class 'unit.tree'. The 'unit.tree' object
#' can then be used to assess the adequacy of phylogenetic models of continuous character
#' evolution. Alternatively, the function can take a 'phylo' object or a 'multiPhylo' object
#' and compute the contrasts on the phylogeny without rescaling. This is only meaningful for
#' downstream analyses if the phylogenies have been rescaled beforehand. Currently,
#' the following object types have been implemented:
#' \itemize{
#'  \item a 'gfit' object returned from fitting a model of continuous character evolution using
#'   \code{\link{geiger::fitContinuous}} in the 'geiger' package.
#'  \item a 'fit.mle' object returned from fitting a model of continuous character evolution
#'   using \code{\link{diversitree::find.mle}} in the 'diversitree' package. As the 'fit.mle' object
#'   does not include all of the information required for creating a 'unit.tree', a second argument
#'   \code{lik} needs to be supplied, providing the likelihood function used in 'find.mle'.
#'  \item a 'gls' object returned from fitting a phylogenetic generalized least squared model
#'   of character correlation using \code{\link{nlme::gls}} in the 'nlme' package.
#'  \item a 'phylo' object. If a 'phylo' object is supplied, the tree is assumed to have been
#'   rescaled previously. A second argument \code{data} must also be provided included the trait
#'   data as a named vector with names equal to the tip.labels of the phylogeny.
#' }
#' 
#'
#' @return a 'unit.tree' object containing (or a list of 'unit.tree' objects,
#' each containing) the following elements:
#' \describe{
#'  \item{phy}{a 'phylo' object. If a model fitted object has been supplied
#' (see details), the 'phylo' object will be rescaled based on fitted model
#' parameters.}
#'  \item{data}{the original comparative data}
#'  \item{pics}{a matrix consisting of the contrasts calculated using the original
#' data on the (rescaled) phylogeny. The matrix has two columns: the "contrasts" and
#' a "variance" for each contrast.}
#' }
#'
#' @export as.unit.tree
#'
#' @seealso \code{\link{ape::pic}}, \code{\link{phy.model.check}}
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
#' ## this creates a 'gfit' object which can be used
#' ## in 'as.unit.tree()'
#' unit.tree.geiger <- as.unit.tree(fit.bm)
#' unit.tree.geiger
#'
#'
#' ## fit Brownian motion model
#' ## using diversitree's find.mle function
#' 
#' ## library(diversitree)
#' ## bmlik <- make.bm(phy, data)
#' ## fit.bm.dt <- find.mle(bmlik, 0.1)
#'
#' ## this creates a 'fit.mle' object which can be used
#' ## in 'as.unit.tree()'
#' ## unit.tree.dt <- as.unit.tree(fit.bm.dt)
#'
#'
#' ## Use pgls to look for a correlation between two traits
#' ## library(nlme)
#' ## t1 <- data
#' ## t2 <- td$data[,"tarsusL"]
#' ## dd <- cbind.data.frame(t1, t2)
#'
#' ## fit gls model with corPagel correlation structure
#' ## fit.gls <- gls(t1~t2, data=dd, correlation=corPagel(phy=phy, value=1))
#'
#' ## this creates a 'gls' object which can be used
#' ## in 'as.unit.tree()'
#' ## unit.tree.gls <- as.unit.tree(fit.gls)
#' 
as.unit.tree <- function(x, ...)
    UseMethod("as.unit.tree")







#' @method as.unit.tree phylo
#' @S3method as.unit.tree phylo
as.unit.tree.phylo <- function(x, data, ...) {
  ## check tree and data to make sure they match
  td <- treedata(phy=x, data=data)
  phy <- td$phy
  data <- td$data

  ## calculate pics
  pics <- pic(data, phy, var.contrasts=TRUE)

  ## append all the object together
  unit.tree <- list(phy=phy, data=data, pics=pics)

  ## change the class of the unit.tree
  class(unit.tree) <- "unit.tree"

  unit.tree
}







#' @method as.unit.tree gfit
#' @S3method as.unit.tree gfit
as.unit.tree.gfit <- function(x, ...) {
  ## use S3 generic modelinfo to pull out the tree, data, parameter
  ## estimates and model type
  obj <- model.info(x, ...)

  ## rescale the phylogeny according to the model
  phy <- make.model.phylo(obj)

  ## build unit.tree from phylo object
  as.unit.tree(phy, obj$data$data)
}








#' @method as.unit.tree gls
#' @S3method as.unit.tree gls
as.unit.tree.gls <- function(x, ...){
    ## use modelinfo to pull out tree data and parameter estimates
    obj <- model.info(x, ...)

    ## rescale the phylogeny according to the model
    phy <- make.model.phylo(obj)

    ## build unit.tree from phylo object
    ## here, data is the residuals
    as.unit.tree(phy, obj$data$data)
}









#' @method as.unit.tree fit.mle
#' @S3method as.unit.tree fit.mle
as.unit.tree.fit.mle <- function(x, ...){
    ## use modelinfo to get tree, data and model parameters
    obj <- model.info(x, ...)

    ## rescale tree according to the model
    phy <- make.model.phylo(obj)

    ## create unit.tree from rescaled phylogeny and data
    as.unit.tree(phy, obj$data$data)
}








#' @title Check if object is a unit.tree
#'
#' @description Utility function for checking if object is of class 'unit.tree'
#'
#' @param x object
#'
#' @return logical
#'
#' @seealso \code{\link{as.unit.tree}}
#'
#' @export is.unit.tree
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
#' is.unit.tree(unit.tree.phy)
#' 
is.unit.tree <- function(x)
    inherits(x, "unit.tree")







## add function for producing error if unit.tree is expected and not provided
assert.is.unit.tree <- function(x){
    if (!inherits(x, "unit.tree"))
        stop(deparse(substitute(x)), " must be a 'unit.tree'")
}
