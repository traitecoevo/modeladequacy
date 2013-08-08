make.sigsqReml <-
function(unit.tree){
	
	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class phylo"))
	}
	
	.sigsqReml <- function(data){
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data
	
	## Take pics
	pics <- pic(data, tree) 
	
	## Get mean of squared pics 
	remlss <- mean(pics^2)	
	
	remlss
	
	}
	
	list(fxn=.sigsqReml, unit.tree=unit.tree)
	
}
