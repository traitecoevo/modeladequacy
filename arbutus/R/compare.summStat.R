## arbutus:::compare.summStat

## fxn for comparing empirical to simulated summary statistics

## calculates p value from quantile

## takes 2 arguments:

## summ.stats.obs -- the observed summary statistics. Can be a vector or a data frame (if samples from posterior). Note that if there are summ stats from multiple unit.trees, it must be the same length as summ.stats.sim

## summ.stats.sim -- the distribution of simulated summary statistics. Input as a data.frame with colnames corresponding to the summary statistics used


## returns a named vector of p.values, with each element corresponding to a statistic used



## NOTE: need to check to make sure I am computing the "p-values" correctly



compare.summStat <- function(summ.stats.obs, summ.stats.sim){
	
	## check to make sure names are the same
	if (!all(names(summ.stats.obs) %in% names(summ.stats.sim)))
		stop("Column names must match between two sets of summary statistics")


	## if multiple columns, reorder to make sure they are in the correct order
        if (ncol(summ.stats.obs) > 1){
            summ.stats.sim <- summ.stats.sim[,names(summ.stats.obs)]
        }
	
	if (nrow(summ.stats.obs) == 1){ # only one set of summ.stats supplied
		
		p.value <- vector()
		
		for (i in 1:ncol(summ.stats.obs)){
			
			## empirical summary stats
			obs <- summ.stats.obs[,i]
			
			## simulated summary stats
			sim <- summ.stats.sim[,i]
			
			## calculate p-value
			p.r <- length(which(sim >= obs))/(length(sim) + 1)
			p.l <- length(which(sim < obs))/(length(sim) + 1)
			p <- min(p.r, p.l)
			
			p.value <- c(p.value, p)
		}
		
		names(p.value) <- colnames(summ.stats.obs)
		
	} else {
	
		## if data.frame of observed summary statistics supplied
		## check to make sure that they are the same size as the
		## simulated summary statistics
		if(nrow(summ.stats.obs) != nrow(summ.stats.sim))
			stop("If inputting multiple observed summary stats, the data.frame must be the same size as that of the simulated summary stats")


                p.value <- vector()

		for (j in 1:ncol(summ.stats.obs)){

                    p.l <- vector()
                    p.r <- vector()

                    ## compute pairwise comparisons
			for (i in 1:nrow(summ.stats.obs)){

                            tmp.l  <- summ.stats.sim[i,j] < summ.stats.obs[i,j]
                            p.l <- c(p.l, tmp.l)
                            tmp.r <- summ.stats.sim[i,j] >= summ.stats.obs[i,j]
                            p.r <- c(p.r, tmp.r)

                        }

                    p.l <- length(which(p.l))/(nrow(summ.stats.obs))
                    p.r <- length(which(p.r))/(nrow(summ.stats.obs))
                    p <- min(p.l, p.r)

                    p.value <- c(p.value, p)
			
		}
		
		names(p.value) <- colnames(summ.stats.obs)
					
	}
	
	## return two-tailed p-values
	pp <- p.value*2

        ## return all info so that summary stats can be plotted
        out <- list(p.values = pp, summ.stats.obs=summ.stats.obs, summ.stats.sim=summ.stats.sim)

        ## define class for print fxn
	class(out) <- c("phy.ss", "list")
        return(out)
			
}




## define print fxn for class phy.ss

print.phy.ss <- function(x){
   cat("Two-tailed p-values for summary statistics: \n")
   cat("\n")
   print(unclass(x$p.values))
}


## function to extract p-values from phy.ss object

pvalueSummStat <- function(x)
    x$p.values



