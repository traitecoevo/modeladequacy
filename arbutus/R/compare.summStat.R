#' @title Compare observed and simulated summary statistics
#'
#' @description Calculates two--tailed p-values from comparing
#' observed to simulated summary statistics in order to evaluate
#' model adequacy for phylogenetic models of continuous character evolution.
#'
#' @param summ.stats.obs a data.frame containing the observed summary statistics
#' with the column names denoting the names of the statistics used.
#'
#' @param summ.stats.sim a data.frame containing the simulated summary statistics
#' with the column names denoting the names of the statistics used.
#'
#' @details The column names must be the same for the observed and simulated
#' summary statistics (though not necessarily in the same order). Both
#' \code{summ.stats.obs} and \code{summ.stats.sim} are designed to be obtained
#' from the function \code{\link{summ.stats}}.
#'
#' If there is only 1 row in \code{summ.stats.obs} (such as the case if ML estimate
#' of model is used), then the function will compare the observed summary statistic(s) to
#' the distribution of simulated simulated summary statistics.
#'
#' If there is more than 1 row in \code{summ.stats.obs} (such as the case if we have
#' a posterior distribution of parameter estimates or the summary statistics are calculated
#' across many trees), then the function requires that the number of rows in
#' \code{summ.stats.obs} and \code{summ.stats.sim} be equal. This is because the function
#' calculates the p-values by performing pairwise comparisons. the function will automatically
#' check that this condition is met and produce an error if it is not.
#'
#' In both cases, the p-values that are returned are two--tailed p-values as
#' the function considers the probability of the observed statistics being larger
#' or smaller than the simulated statistics.
#'
#' The distribution of simulated summary statistics along with the observed summary
#' statistic using a generic \code{\link{plot}} function (see \code{\link{plot.phy.ss}})
#' for more details. For the purposes of plotting, the observed and simulated summary
#' statistics used as arguments are returned with the p-values.
#'
#' The p-values can be extracted from the output of \code{compare.summ.stats} using the
#' function \code{\link{pval.summ.stats}}
#'
#' @return an object of class 'phy.ss' which is a list containing the following objects:
#' \itemize{
#'  \item{\code{p.values}: }{Two--tailed p-values for all summary statistics}
#'  \item{\code{summ.stats.obs}: }{The input data frame of observed summary statistics}
#'  \item{\code{summ.stats.sim}: }{The input data frame of simulated summary statistics}
#' }
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{plot.phy.ss}}, \code{\link{pval.summ.stats}},
#' \code{\link{phy.model.check}}
#'
#' @export compare.summ.stats
#'
#' @author Matt Pennell
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## calculate default summary stats on observed data
#' ss.obs <- summ.stats(unit.tree, stats=NULL)
#'
#' ## simulate data on unit.tree
#' sims <- sim.char.unit(unit.tree, nsim=10)
#'
#' ## calculate default summary stats on simulated data
#' ss.sim <- summ.stats(sims, stats=NULL)
#'
#' ## compare simulated to observed summary statistics
#' res <- compare.summ.stats(ss.obs, ss.sim)
#'
#' ## extract p-values
#' pval.summ.stats(res)
#'
#' ## plot results
#' plot(res)
#' 
compare.summ.stats <- function(summ.stats.obs, summ.stats.sim){
	
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










#' @title Print summary statistics
#'
#' @description Generic function for printing 'phy.ss' objects
#' from \code{\link{compare.summ.stats}}
#'
#' @param x a 'phy.ss' object
#'
#' @return prints p.values to screen
#'
#' @author Matt Pennell
#'
#' @method print phy.ss
#' @S3method print phy.ss
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#'
#' ## fit BM model using fitContinuous
#' f <- fitContinuous(phy, dat, model="BM", control=list(niter=10))
#'
#' ## use phy.model.check to assess model adequacy
#' p <- phy.model.check(f, nsim=10)
#'
#' ## print output
#' p
print.phy.ss <- function(x){
   cat("Two-tailed p-values for summary statistics: \n")
   cat("\n")
   print(unclass(x$p.values))
}






#' @title Extract p--values for summary statistics
#'
#' @description Utility function for extracting p-values from the output
#' of \code{\link{compare.summ.stats}}
#'
#' @param x a 'phy.ss' object from the function \code{\link{compare.summ.stats}}
#'
#' @return a named vector of two-tailed p-values
#'
#' @export pval.summ.stats
#'
#' @seealso \code{\link{compare.summ.stats}}
#'
#' @author Matt Pennell
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## calculate default summary stats on observed data
#' ss.obs <- summ.stats(unit.tree, stats=NULL)
#'
#' ## simulate data on unit.tree
#' sims <- sim.char.unit(unit.tree, nsim=10)
#'
#' ## calculate default summary stats on simulated data
#' ss.sim <- summ.stats(sims, stats=NULL)
#'
#' ## compare simulated to observed summary statistics
#' res <- compare.summ.stats(ss.obs, ss.sim)
#'
#' ## extract p-values
#' pval.summ.stats(res)
#'
pval.summ.stats <- function(x)
    x$p.values



