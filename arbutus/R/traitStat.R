## arbutus:::traitStat

## Fxn for calculating summary statistics across a unit.tree or a set of unit.trees

## 3 arguments
## unit.tree -- can be either a single unit.tree or a list of unit.trees

## data -- default is NULL, meaning that the data from the unit.tree(s) is used. If data != NULL, it will apply the same dataset to all unit.trees

## summ.stats -- a named list of summary statistics. Default is NULL, meaning that the fxn uses the set of 6 default summary statistics


## returns a named data.frame with all of the summary statistics


traitStat <- function(unit.tree, data=NULL, stats=NULL){
	if (is.null(stats))
		stats <- defSummStats()
		
	res <- data.frame()
	
	## compute all stats across all unit.trees 
	for (i in 1:length(unit.tree)){
		tmp <- lapply(stats, function(x) x(unit.tree[[i]], data))
		tmp2 <- unlist(tmp)
		res <- rbind(res,tmp2)
	}
	
	res 	
	
}