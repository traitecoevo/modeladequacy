## quercus:::sigsqReml

## A function for computing the REML estimate of sigma^2
## Equal to the mean of the squared contrasts
## Takes a unit.tree object
## If data == NULL, uses the data supplied by unit.tree
## If data != NULL, uses input data

sigsqReml <- function(unit.tree, data=NULL){
	
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)
		
	mean(unit.tree$pics[,"contrasts"]^2)	
		
}