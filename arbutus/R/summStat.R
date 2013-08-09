## arbutus:::summStat

## Fxn for calculating summary statistics across a unit.tree or a set of unit.trees

## 2 arguments
## unit.tree -- can be either a single unit.tree or a list of unit.trees

## summ.stats -- a named list of summary statistics. Default is NULL, meaning that the fxn uses the set of 6 default summary statistics


## returns a named data.frame with all of the summary statistics

## Note: data from unit.trees will be used


summStat <- function(unit.tree, stats=NULL){
	if (is.null(stats))
		stats <- defSummStats()
	
	## single unit.tree
	if (inherits(unit.tree, "unit.tree")){
		res <- lapply(stats, function(x) x(unit.tree))
		res <- do.call(cbind.data.frame, res)
		
	} else { ## multiple unit.trees	
				
		res <- data.frame()
	
		## compute all stats across all unit.trees 
		for (i in 1:length(unit.tree)){
			tmp <- lapply(stats, function(x) x(unit.tree[[i]]))
			tmp2 <- do.call(cbind, tmp)
			res <- rbind(res,tmp2)
			
		}
	}	
	
	res 	
	
}