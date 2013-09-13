as.unit.tree <- function(x, ...)
  UseMethod("as.unit.tree")

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

as.unit.tree.gfit <- function(x, ...) {
  info <- modelinfo(x)
  as.unit.tree(rescale(info), info$data$data)
}

as.unit.tree.fit.mle <- function(x, lik, ...) {
  info <- modelinfo(x, lik)
  as.unit.tree(rescale(info), info$data$data)
}

## This is temporary, and will depend on Matt's changes to the
## rescaling functions.  Only handles BM for now.
rescale <- function(obj) {
  if (obj$type == "BM") {
    phy <- obj$data$phy
    sigsq <- obj$pars[["sigsq"]]
    se <- obj$pars[["SE"]]

    phy$edge.length <- phy$edge.length * sigsq
    if (!is.null(se)) {
      tips <- phy$edge[,2] <= Ntip(phy)
      phy$edge.length[tips] <- phy$edge.length[tips] + se^2
    }

    phy
  } else {
    stop("Not yet implemented")
  }
}

is.unit.tree <- function(x)
  inherits(x, "unit.tree")
