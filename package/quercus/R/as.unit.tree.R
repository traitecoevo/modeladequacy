## quercus:::as.unit.tree()

## takes a phylo object (which has been previously rescaled) and a data frame
## calculates the contrasts and variance
## appends both the data frame and the contrasts to the phylo object
## creates a new class of object called a 'unit.tree'

as.unit.tree <- function(x, data){
	
	## check to make sure object is of class 'phylo'
	if (!inherits(x, "phylo"))
		stop("rescaled tree must be of class `phylo'")
	
	## check tree and data to make sure they match	
	td <- treedata(phy=x, data=data)
	phy <- td$phy
	data <- td$data
	
	## calculate pics
	pics <-	pic(data, phy, var.contrasts=TRUE)
	
	## append all the object together
	unit.tree <- list(phy=phy, data=data, pics=pics)
	
	## change the class of the unit.tree
	class(unit.tree) <- c("unit.tree", "phylo")
	
	unit.tree

}	
		
