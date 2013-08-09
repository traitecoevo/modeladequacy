## arbutus:::nodeAge

## Fxn for computing the branching times of a phylogeny
## Similar to ape's branching.times() fxn but works with non-ultrametric trees
## Used internally to compute the summary statistic in which
## a linear model is fit between the node height and contrasts

nodeAge <- function(x}{
	
	if (!inherits(x, "phylo"))
		stop("object must be of class 'phylo'")
		
	n <- length(x$tip.label)
	N <- dim(phy$edge)[1]
	
	nodes <- numeric(x$Nnode)
	
	int <- which(x$edge[,2] > n)
	
	for (i in int){
		nodes
	}
		
	
}