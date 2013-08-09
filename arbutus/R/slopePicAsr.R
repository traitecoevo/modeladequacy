## arbutus:::slopePicAsr

## Fxn for computing the slope between the contrasts and the inferred ancestral state at the node
## Used as a measure of variation with respect to trait value

## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data



## NOTES:
##
## 1. Use "pic" method in ape. Not sure if this is appropriate or not.
 

slopePicAsr <- function(unit.tree, data=NULL){
	
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	## reconstruct ancestors using pic method
	asr <- ace(unit.tree$data, unit.tree$phy, method="pic")$ace
	
	## absolute value of pics
	abs.pic <- abs(unit.tree$pics[,"contrasts"])
	
	## fit linear model
	ca <- lm(abs.pic ~ asr)
	
	slope <- ca$coefficients["asr"]
	
	as.numeric(slope)
		
}
