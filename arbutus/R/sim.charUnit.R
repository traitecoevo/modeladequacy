#' @title Simulate character evolution on unit tree
#'
#' @description Simulates continuous characters along a unit tree according to a Brownian motion
#'   process with rate 1.
#'
#' @param unit.tree a 'unit.tree' object
#' @param nsim the number of datasets to simulate (if single unit.tree is provided; see below)
#' 
#' @return A list of unit.trees. In each unit.tree, the phylogeny will be the same as that of the input unit.tree(s)
#'   but the data and the contrasts will be replaced by the simulated datasets.
#' 
#' @details If a single unit.tree is supplied, the function will simulate \code{nsim} datasets under a Brownian motion
#'   process with a rate of 1. If a lits of unit.trees supplied (such as those derived from a Bayesian analysis)
#'   \code{sim.char.unit} will simulate a single data set on each tree in the list and the \code{nsim} argument will
#'   be ignored.
#' 
#' @export sim.char.unit
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#'
#' ## build unit.tree object
#' unit.tree <- as.unit.tree(phy, data=dat)
#'
#' ## simulate 2 datasets on unit tree
#' sims <- sim.char.unit(unit.tree, nsim=2)
#' 
#' sims
#' 
sim.char.unit <- function(unit.tree, nsim=1000) {
  if (inherits(unit.tree, "unit.tree")) {
    phy <- unit.tree$phy
    dat <- sim.char.std.bm(phy, nsim)
    lapply(seq_len(nsim), function(i) as.unit.tree(phy, dat[,i]))
  } else if (is.list(unit.tree)) {
    ## Simulate one data set per tree
    lapply(unit.tree, function(u) sim.char.unit(u, nsim=1)[[1]])
  } else {
    stop("'unit.tree' must be a single unit tree or a list of unit.trees")
  }
}

## Actual character simulation code imported from diversitree.  This
## is optimised for the case where we have BM with rate 1.
sim.char.std.bm <- function(tree, nsim=1, x0=0) {
  edge <- tree$edge
  idx <- seq_len(max(edge))
  n.tip <- length(tree$tip.label)
  root <- n.tip + 1
  is.tip <- idx <= n.tip
  children <- get.children(edge, n.tip)
  order <- rev(get.ordering(children, is.tip, root))
  len <- tree$edge.length[match(idx, edge[, 2])]

  len[root] <- 0 # set root length to zero to avoid warning.

  n.edge <- length(len)
  dy <- matrix(rnorm(n.edge * nsim, 0, rep(sqrt(len), nsim)),
               n.edge, nsim)

  y <- matrix(NA, n.edge, nsim)
  y[root,] <- x0

  for (i in order) {
    j <- children[[i]]
    y[j,] <- y[rep.int(i, length(j)),] + dy[j,]
  }

  y.tip <- y[is.tip,,drop=FALSE]
  rownames(y.tip) <- tree$tip.label

  y.tip
}

## These two could probably be replaced by using ape's reorder, but
## are also taken from diversitree.
get.ordering <- function(children, is.tip, root) {
  todo <- list(root)
  i <- root
  repeat {
    kids <- unlist(children[i])
    i <- kids[!is.tip[kids]]
    if (length(i) > 0)
      todo <- c(todo, list(i))
    else
      break
  }
  as.vector(unlist(rev(todo)))
}

## Note that unlike diversitree's version, this returns a list and
## allows multifurcations.
get.children <- function(edge, n.tip) {
  ## To construct the children vector, this is what I used to do:
  ##   lapply(idx[!is.tip], function(x) edge[edge[,1] == x,2])
  ## But this is slow for large trees.  This is faster:
  ## children <- split(edge[,2], edge[,1])
  ## Surprisingly, most of the time is in coercing edge[,1] into a
  ## factor.
  x <- as.integer(edge[,1])
  levels <- as.integer((n.tip+1):max(edge[,1]))
  f <- match(x, levels)
  levels(f) <- as.character(levels)
  class(f) <- "factor"
  children <- split(edge[,2], f)
  names(children) <- NULL

  c(rep(list(integer(0)), n.tip), children)
}
