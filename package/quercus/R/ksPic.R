## quercus:::ksPic

## A function for comparing the distribution of the contrasts
## to a normal distribution with mean 0 and sd = sqroot(mean(pics^2))
## uses D-statistic from KS test

ksPic <- function(unit.tree, data=NULL){
	
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)
		
	## simulate a null distribution
	sd <- sqrt(mean(unit.tree$pics[,"contrasts"]^2))	
	nd <- rnorm(10000, mean=0, sd=sd)
	
	## KS test
	ksbm <- ks.test(pics[,"contrasts"], nd)$statistic
	
	as.numeric(ksbm)	

}