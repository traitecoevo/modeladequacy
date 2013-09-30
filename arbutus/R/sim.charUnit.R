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
#' @seealso \code{\link{geiger::sim.char}} which this function wraps
#' 
#' @export sim.char.unit
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#'
#' ## build unit.tree object
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## simulate 2 datasets on unit tree
#' sims <- sim.char.unit(unit.tree, nsim=2)
#' 
#' sims
#' 
sim.char.unit <- function(unit.tree, nsim=1000){
	if (inherits(unit.tree, "unit.tree")){ ## simulate on one tree
            phy <- unit.tree$phy
            dat <- sim.char(phy, par=1, nsim=nsim, model="BM")[,,]
            ut <- lapply(1:nsim, function(x) as.unit.tree(phy, dat[,x]))
                
            } else { ## simulate one simulation per tree in list
             ## first check to make sure first element is a unit.tree
            if (!is.unit.tree(unit.tree[[1]]))
                stop("unit.tree must be either a single unit tree or a list of unit.trees")

                multi.phy <- lapply(unit.tree, function(x) return(x$phy))
                multi.dat <- lapply(multi.phy, function(x) return(sim.char(x, par=1, model="BM")[,,]))
                ut <- lapply(1:length(multi.phy), function(x) as.unit.tree(multi.phy[[x]], multi.dat[[x]]))
            }
                                    
	ut
	
}
