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
#' The p-values can be extracted from the output of \code{\link{compare.summ.stats}} using the
#' function \code{\link{pval.summ.stats}}
#'
#' @return an object of class \code{phy.ss} which is a list containing the following objects:
#' \describe{
#'  \item{p.values}{Two--tailed p-values for all summary statistics}
#'  \item{summ.stats.obs}{The input data frame of observed summary statistics}
#'  \item{summ.stats.sim}{The input data frame of simulated summary statistics}
#' }
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{plot.phy.ss}}, \code{\link{pval.summ.stats}},
#' \code{\link{phy.model.check}}
#'
#' @export compare.summ.stats
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, data=dat)
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
#' @param x a \code{phy.ss} object
#' @param ... additional arguments to be passed to print
#'
#' @return prints p.values to screen
#'
#'
#' @method print phy.ss
#' @S3method print phy.ss
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#'
#' ## use phy.model.check to assess model adequacy
#' ## pretend as if tree has already been rescaled
#' p <- phy.model.check(phy, data=dat, nsim=10)
#'
#' ## print output
#' p
#' 
print.phy.ss <- function(x, ...){
   cat("Two-tailed p-values for summary statistics: \n")
   cat("\n")
   print(unclass(x$p.values))
}






#' @title Extract p--values for summary statistics
#'
#' @description Utility function for extracting p-values from the output
#' of \code{\link{compare.summ.stats}}
#'
#' @param x a \code{phy.ss} object from the function \code{\link{compare.summ.stats}}
#'
#' @return a named vector of two-tailed p-values
#'
#' @export pval.summ.stats
#'
#' @seealso \code{\link{compare.summ.stats}}
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, data=dat)
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





#' @title Multivariate measure of model adequacy
#'
#' @description Computes Mahalanobis distance between the observed summary statistics and
#' the simulated summary statistics as a multivariate measure of model fit
#'
#' @param x a \code{phy.ss} object inherited from \code{\link{compare.summ.stats}}
#'
#' @details This function computes the Mahalanobis distance between the observed and simulated
#' summary statistics. The Mahalanobis distance (see \code{\link[stats]{mahalanobis}} is a unit-less and
#' scale-invariant of the distance
#' between a single data point (our observed summary statistics) and a common point (here the mean of the
#' simulated summary statistics), taking into account the covariance between the summary statistics from
#' the simulated data. It assumes that the distribution of summary statistics is multivariate normal.
#' For the default summary statistics (see \code{\link{def.summ.stats}}, this condition should be met -- the
#' exception being the \code{\link{ks.contrast}} statistic, which as it is bounded at 0, will not be. As a result,
#' if \code{\link{ks.contrast}} is included in the set of summary statistic, this function will take the
#' natural log of the values before computing the Mahalanobis distance. All other summary statistics will
#' be assumed to be normally distributed and used as is.
#'
#' While the Mahalanobis distance may be a useful summary measure in some circumstances, we recommend checking
#' the summary statistics individually for a number of reasons. First, our procedure for calculating the p-values
#' for the summary statistics is general and does not depend on assumptions regarding the distribution of values.
#' Second, the interpretation from the individual p-values is much more clear from the perspective of either
#' posterior predictive or parametric bootstrapping theory. Third, and most importantly, the fact that some summary
#' statistics capture the variation in the data better than others provides useful information as to how and why the
#' model is inadequate.
#'
#' If only one set of observed
#' summary statistics are available (e.g. from fitting a model using maximum likelihood to a single tree), a single distance
#' will be returned. If multiple sets of summary statistics are available (e.g. from fitting a model using a Bayesian MCMC),
#' the function will return a distribution of distances.
#'
#' @return the Mahalanobis distance between the observed and simulated summary statistic. 
#'
#' @export mv.summ.stats
#'
#' @seealso \code{\link{compare.summ.stats}}, \code{\link[stats]{mahalanobis}}
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, data=dat)
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
#' ## calculate Mahalanobis distance
#' mv.summ.stats(res)
#'
mv.summ.stats <- function(x){
    if (!inherits(x, "phy.ss"))
        stop("x must be a 'phy.ss' object")
    ## log ks d-statistic
    if ("ks.contrast" %in% names(pval.summ.stats(x))){
        x$summ.stats.obs[,"ks.contrast"] <- log(x$summ.stats.obs[,"ks.contrast"])
        x$summ.stats.sim[,"ks.contrast"] <- log(x$summ.stats.sim[,"ks.contrast"])
    }
    obs <- as.matrix(x$summ.stats.obs)
    sim <- as.matrix(x$summ.stats.sim)
    cv <- cov(sim)

    m <- mahalanobis(x=obs, center=colMeans(sim), cov=cv)
    m
}
