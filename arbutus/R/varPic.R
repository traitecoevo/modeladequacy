## arbutus:::varPic

## A fxn for computing the variance in pic estimates
## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data

varPic <- function(unit.tree, data=NULL) {
	
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)
		
  	var(abs(unit.tree$pics[,"contrasts"]))

}
