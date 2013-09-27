## arbutus:::sim.char.unit

## fxn for simulating BM with rate 1 along a unit tree
## repeats n times

## each simulation generates a unit.tree with the simulated data

## 2 arguments:
## unit.tree -- the unit.tree across which to simulate (can just be of class 'phylo' as data or pics are not used)

## nsim -- the number of simulations

## Note that if it is given one unit.tree, the fxn simulates nsim datasets on the unit tree
## If given a list of unit.trees, the fxn simulates one dataset per tree and ignores the nsim command


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
