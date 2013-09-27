## arbutus:::as.unit.tree

## x can be either a fitted object or else a phylo object
as.unit.tree <- function(x, ...)
    UseMethod("as.unit.tree")

## x is of class phylo
## data is also required
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



## x is of class gfit (geiger fitted object)
as.unit.tree.gfit <- function(x, ...) {
  ## use S3 generic modelinfo to pull out the tree, data, parameter
  ## estimates and model type
  obj <- model.info(x, ...)

  ## rescale the phylogeny according to the model
  phy <- make.model.phylo(obj)

  ## build unit.tree from phylo object
  as.unit.tree(phy, obj$data$data)
}

## x is of class gls (gls fitted model)
as.unit.tree.gls <- function(x, ...){
    ## use modelinfo to pull out tree data and parameter estimates
    obj <- model.info(x, ...)

    ## rescale the phylogeny according to the model
    phy <- make.model.phylo(obj)

    ## build unit.tree from phylo object
    ## here, data is the residuals
    as.unit.tree(phy, obj$data$data)
}

## x is of class fit.mle (diversitree fitted object)
as.unit.tree.fit.mle <- function(x, ...){
    ## use modelinfo to get tree, data and model parameters
    obj <- model.info(x, ...)

    ## rescale tree according to the model
    phy <- make.model.phylo(obj)

    ## create unit.tree from rescaled phylogeny and data
    as.unit.tree(phy, obj$data$data)
}

## utility fxn to check if object is of class unit.tree
is.unit.tree <- function(x)
    inherits(x, "unit.tree")

## add function for producing error if unit.tree is expected and not provided
assert.is.unit.tree <- function(x){
    if (!inherits(x, "unit.tree"))
        stop(deparse(substitute(x)), " must be a 'unit.tree'")
}
