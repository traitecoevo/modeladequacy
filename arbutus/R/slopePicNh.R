## arbutus:::slopePicNh

## Fxn for computing the slope between the contrasts and the height of the node at which the contrast is computed
## Used as a measure of variation with respect to time

## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data

# NOTES:
##
## 1. Need to swap out branching.times() with a function that works for non-ultrametric trees (working on the fxn nodeAge to do this)

slopePicNh <- function(unit.tree, data=NULL){
	
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)
	
	## node heights (see note 1)
	nh <- branching.times(unit.tree$phy)
	
	## absolute value of the pics	
	abs.pic <- abs(unit.tree$pics[,"contrasts"])
	
	## fit a linear model
	
	ch <- lm(abs.pic ~ nh)
	
	slope <- ch$coefficients["nh"]
	
	as.numeric(slope)	

}