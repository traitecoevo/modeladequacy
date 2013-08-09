## arbutus:::simcharUnit

## fxn for simulating BM with rate 1 along a unit tree
## repeats n times

## each simulation generates a unit.tree with the simulated data

## 2 arguments:
## unit.tree -- the unit.tree across which to simulate (can just be of class 'phylo' as data or pics are not used)

## nsim -- the number of simulations


sim.charUnit <- function(unit.tree, nsim){
	if (!inherits(unit.tree, "unit.tree"))
		stop("unit.tree must be of class 'unit.tree'")
		
	phy <- unit.tree$phy	
	
	dat <- sim.char(phy, par=1, nsim=nsim, model="BM")[,,]
	
	ut <- lapply(1:nsim, function(x) as.unit.tree(phy, dat[,x]))
	
	ut
	
}