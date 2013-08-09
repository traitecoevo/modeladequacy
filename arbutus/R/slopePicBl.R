## arbutus:::slopePicBl

## Fxn for computing the slope between the contrasts and their standard deviations
## Used as a measure of variation with respect to branch length

## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data

slopePicBl <- function(unit.tree, data=NULL){
	
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	## sd of the pics
	sd.pic <- sqrt(unit.tree$pics[,"variance"])
	
	## absolute value of the pics
	abs.pic <- abs(unit.tree$pics[,"contrasts"])
	
	## fit a linear model
	cv <- lm(abs.pic ~ sd.pic)
	
	as.numeric(cv$coefficients["sd.pic"])
	
}