## plotting functions for arbutus summary statistics


#' @title Plot contrasts versus their variances
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' their expected variances.
#'
#' @param x a model fitted object or \code{phylo} object that can be passed to \code{\link{as.unit.tree}}.
#'
#' @param col colours used for plotting
#'
#' @param ... additional arguments to be passed to \code{\link{as.unit.tree}}
#'
#' @details This function generates a plot of the observed contrasts versus their
#' expected variances (see \code{\link[ape]{pic}}). The slope of this plot is equal to
#' the summary statistic \code{\link{cor.contrast.var}}
#' and is used as a default summary statistic in \code{\link{summ.stats}}.
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single dataset (does not integrate across multiple trees).
#'
#' \code{contrast.var.plot} uses the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @seealso \code{\link{as.unit.tree}}, \code{\link{summ.stats}}, \code{\link{cor.contrast.var}}
#'
#' @export contrast.var.plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' \dontrun{
#' require(geiger)
#' ## fit Brownian motion model
#' ## using geiger's fitContinuous function
#' fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
#'                                  control=list(niter=10))
#'
#' ## plot contrasts versus their variances
#' contrast.var.plot(fit.bm)
#' }
#' 
contrast.var.plot <- function(x, col=c("dodgerblue4", "darkblue"), ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## create data frame for pics
    dat <- as.data.frame(ut$pics)
    dat$contrasts <- abs(dat$contrasts)

    .e <- environment()
    p <- ggplot(dat, aes(variance, contrasts), environment=.e)
    p <- p + stat_smooth(method="lm", fill=col[1], colour=col[2], size=1, alpha=0.1)
    p <- p + geom_point()
    p <- p + theme_bw()
    p <- p + xlab("E[Var] of contrast")
    p <- p + ylab("Contrast")
    p <- p + ggtitle("Variation relative to branch lengths")

    print(p)
}














#' @title Plot contrasts versus their inferred ancestral states
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' their inferred ancestral state.
#'
#' @param x a model fitted object or \code{phylo} object that can be passed to \code{\link{as.unit.tree}}.
#'
#' @param col colours used for plotting
#'
#' @param ... additional arguments to be passed to \code{\link{as.unit.tree}}
#'
#' @details This function generates a plot of the observed contrasts versus their
#' inferred ancestral state (see \code{\link[ape]{pic}}, \code{\link[ape]{ace}}).
#' Note that the ancestral states are estimated using the 'pic' method, which
#' differs somewhat from the ML estimate of the ancestral states.
#' The slope of this plot is equal to
#' the summary statistic \code{\link{cor.contrast.asr}}
#' and is used as a default summary statistic in \code{\link{summ.stats}}. 
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single dataset (does not integrate across multiple trees).
#'
#' \code{contrast.asr.plot} uses the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @seealso \code{\link{as.unit.tree}}, \code{\link{summ.stats}}, \code{\link{cor.contrast.asr}}
#'
#' @export contrast.asr.plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' \dontrun{
#' require(geiger)
#' ## fit Brownian motion model
#' ## using geiger's fitContinuous function
#' fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
#'                                  control=list(niter=10))
#'
#' ## plot contrasts versus their variances
#' contrast.asr.plot(fit.bm)
#' }
#' 
contrast.asr.plot <- function(x, col=c("dodgerblue4", "darkblue"), ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## extract pics
    pics <- abs(ut$pics[,"contrasts"])

    ## get tree and data, calculate ancestral states
    asr <- ace(ut$data, ut$phy, method="pic")$ace

    ## create a dataframe
    dat <- cbind.data.frame(asr, pics)

    .e <- environment()
    p <- ggplot(dat, aes(asr, pics), environment=.e)
    p <- p + stat_smooth(method="lm", fill=col[1], colour=col[2], size=1, alpha=0.1)
    p <- p + geom_point()
    p <- p + theme_bw()
    p <- p + xlab("Ancestral state")
    p <- p + ylab("Contrast")
    p <- p + ggtitle("Variation relative to trait value")

    print(p)
}












#' @title Plot contrasts versus their node height
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' the node height (i.e. depth) at which they are computed
#'
#' @param x a model fitted object or \code{phylo} object that can be passed to \code{\link{as.unit.tree}}.
#'
#' @param col colours used for plotting
#'
#' @param ... additional arguments to be passed to \code{\link{as.unit.tree}}
#'
#' @details This function generates a plot of the observed contrasts versus their
#' the node height at which they are computed (see \code{\link[ape]{pic}}).
#' The slope of this plot is equal to
#' the summary statistic \code{\link{cor.contrast.nh}}
#' and is used as a default summary statistic in \code{\link{summ.stats}}. 
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single dataset (does not integrate across multiple trees).
#'
#' \code{contrast.nh.plot} uses the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @seealso \code{\link{as.unit.tree}}, \code{\link{summ.stats}}, \code{\link{cor.contrast.nh}}
#'
#' @export contrast.nh.plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' \dontrun{
#' require(geiger)
#' ## fit Brownian motion model
#' ## using geiger's fitContinuous function
#' fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
#'                                  control=list(niter=10))
#'
#' ## plot contrasts versus their variances
#' contrast.nh.plot(fit.bm)
#' }
#' 
contrast.nh.plot <- function(x, col=c("dodgerblue4", "darkblue"), ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## extract pics
    pics <- abs(ut$pics[,"contrasts"])

    ht <- edge.height(ut$phy)
    N <- Ntip(ut$phy)
    nh <- ht$end[c((N + 1):nrow(ht))]
    names(nh) <- rownames(ht[c((N + 1):nrow(ht)),])

    ## create a dataframe
    dat <- cbind.data.frame(nh, pics)

    .e <- environment()
    p <- ggplot(dat, aes(nh, pics), environment=.e)
    p <- p + stat_smooth(method="lm", fill=col[1], colour=col[2], size=1, alpha=0.1)
    p <- p + geom_point()
    p <- p + theme_bw()
    p <- p + xlab("Node height")
    p <- p + ylab("Contrast")
    p <- p + ggtitle("Variation relative to time")

    print(p)
}




#' @title Plot summary statistics for model adequacy
#'
#' @description Plots the distribution of simulated summary statistics
#' along with the observed summary statistics (either a point estimate or
#' a distribution of values) for visual inspection of model adequacy
#'
#' @param x a 'phy.ss' object inherited from \code{\link{phy.model.check}}
#' (or, alternatively \code{\link{compare.summ.stats}})
#'
#' @param colour a vector of two colours to be used in the plot; the first is for the observed summary
#' statistics and the second for the simulated summary statistics.
#' 
#' @param ... additional arguments to be passed to \code{\link{plot}}. This is currently non-functional. 
#'
#' @details This function is designed to give researchers a quick visual inspection of
#' model adequacy by plotting the simulated distribution of summary statistics along
#' with the empirical estimate(s). If only a point estimated of the observed summary statistics
#' is input (such as from estimating the parameters on a single tree using maximum likelihood),
#' the function will plot the simulated summary statistics with a vertical line representing the
#' observed summary statistic. If a distribution of observed summary statistics is used, the function
#' will plot both distributions (simulated and observed) on the same plot. This function is not
#' designed to be flexible. If researchers want to create their own plots in a different style, all
#' the necessary information is available in the 'phy.ss' object output by both \code{\link{phy.model.check}}
#' and \code{\link{compare.summ.stats}} and can be indexed by \code{x$summ.stats.obs} and
#' \code{x$summ.stats.sim}.
#'
#' This plotting function using the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @method plot phy.ss
#' @S3method plot phy.ss
#'
#' @seealso \code{\link{phy.model.check}}, \code{\link{compare.summ.stats}}
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' 
#' p <- phy.model.check(phy, data=dat, nsim=20)
#'
#' plot(p)
#' 
plot.phy.ss <- function(x, colour=NULL, ...){

    ## check to see if colours were supplied
    ## assigned defaults if not
    if (is.null(colour))
        colour <- c("darkblue", "dodgerblue4")
    
    ## require grid and gridExtra
    require(grid)
    require(gridExtra)
    ## get observed summary stats
    ss.obs <- x$summ.stats.obs
    ss.sim <- x$summ.stats.sim

    ## check if ss.obs is a single value or a distribution
    if (nrow(ss.obs) == 1){ ## historgram plot
        stats <- names(ss.obs)

        p <- lapply(stats, function(x) phy.ss.singleplot(ss.obs[x], ss.sim, colour=colour))

        do.call(grid.arrange, p)
    } else {
        stats <- names(ss.obs)

        p <- lapply(stats, function(x) phy.ss.distplot(ss.obs, ss.sim, x, colour=colour))

        ## dummy plot to get legend with correct colors and style
        legend <- phy.ss.legend(phy.ss.dummyplot(colour))
        lwidth <- sum(legend$width)

        res <- do.call(arrangeGrob, p)
        grid.newpage()
        grid.arrange(res, legend, widths=unit.c(unit(1, "npc") - lwidth, lwidth), ncol=2)
    }
}


                        
## for a single summary statistic                   
phy.ss.singleplot <- function(ss, ss.sim, colour){

    stat <- names(ss)
    range <- range(ss.sim[,stat])
    bins <- (range[2] - range[1]) * 40/length(ss.sim[,stat])

    .e <- environment()
    p <- ggplot(ss.sim, aes(x=ss.sim[,stat]), environment =.e)
    p <- p + geom_histogram(binwidth=bins, alpha=0.6, fill=colour[2], aes(y=..density..))
    p <- p + xlab(stat)
    p <- p + theme_bw()
    p <- p + geom_vline(x=as.numeric(ss), colour=colour[1], size=1, alpha=0.6)

    return(p)
}




## for a distribution of summary statistics
phy.ss.distplot <- function(ss.obs, ss.sim, stat, colour){


    ## bind the two distributions together
    res <- c(ss.obs[,stat], ss.sim[,stat])
    idx <- c(rep("Empirical", length(ss.obs[,stat])), rep("Simulated", length(ss.sim[,stat])))
    dat <- cbind.data.frame(res, idx)
    colnames(dat) <- c("summ.stats", "type")

    range <- range(dat[,"summ.stats"])
    range <- range[2] - range[1]
    bins <- range * 40/length(dat[,"summ.stats"])

    .e <- environment()
    p <- ggplot(dat, aes(x=summ.stats, fill=type), environment = .e)
    p <- p + geom_histogram(binwidth=bins, position="identity", alpha=0.6, aes(y=..density..))
    p <- p + scale_fill_manual(values=colour)
    p <- p + xlab(stat)
    p <- p + theme_bw()
    p <- p + theme(legend.position="none")

    return(p)
}




## internal fxn used to get the legend
## modified from code written by Hadley W.

phy.ss.legend <- function(plot){
    tmp <- ggplot_gtable(ggplot_build(plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

## dummy fxn to quickly get legend
phy.ss.dummyplot <- function(colour){
    x <- rnorm(10)
    y <- c(rep("Empirical", 5), rep("Simulated", 5))
    dat <- cbind.data.frame(x,y)
    colnames(dat) <- c("summ.stats", "type")

    .e <- environment()
    p <- ggplot(dat, aes(x=summ.stats, fill=type), environment = .e)
    p <- p + geom_histogram(binwidth=0.1, alpha =0.6, position="identity", aes(y=..density..))
    p <- p + scale_fill_manual(values=colour)
    p <- p + xlab("dummy")
    p <- p + theme_bw()
    p <- p + theme(legend.title=element_blank())

    return(p)
}
    




                       
    
    
        

        


  
