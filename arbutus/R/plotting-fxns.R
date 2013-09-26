## plotting functions for arbutus summary statistics



## plot slope of contrasts versus variance
## not really that flexible at the moment
## only work for single unit.trees (not posterior distributions)


## takes model fitted objects which will be rerouted through as.unit.tree

contrastVar.plot <- function(x, col=c("dodgerblue4", "darkblue"), ...){

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




## plot slope of the contrasts versus their inferred ancestral state
contrastAsr.plot <- function(x, col=c("dodgerblue4", "darkblue"), ...){

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




## plot of contrasts versus node height

contrastNh.plot <- function(x, col=c("dodgerblue4", "darkblue"), ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## extract pics
    pics <- abs(ut$pics[,"contrasts"])

    ## get tree and get branching times
    nh <- branching.times(ut$phy)

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




## plotting function for summary stats
## compare observed to simulated

## takes a phy.ss object

plot.phy.ss <- function(x, colour=c("dodgerblue4", "darkblue")){

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
    p <- p + geom_histogram(binwidth=bins, alpha=0.5, fill=colour[1], aes(y=..density..))
    p <- p + xlab(stat)
    p <- p + theme_bw()
    p <- p + geom_vline(x=as.numeric(ss), colour=colour[2], size=1, alpha=0.5)

    return(p)
}




## for a distribution of summary statistics
phy.ss.distplot <- function(ss.obs, ss.sim, stat, colour){


    ## bind the two distributions together
    res <- c(ss.obs[,stat], ss.sim[,stat])
    idx <- c(rep("Empirical", length(ss.obs[,stat])), rep("Simulated", length(ss.sim[,stat])))
    dat <- cbind.data.frame(res, idx)
    colnames(dat) <- c("summ.stats", "type")

    .e <- environment()
    p <- ggplot(dat, aes(x=summ.stats), environment = .e)
    p <- p + geom_density(alpha =0.5, aes(fill=type))
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
    p <- ggplot(dat, aes(x=summ.stats), environment = .e)
    p <- p + geom_density(alpha =0.5, aes(fill=type))
    p <- p + scale_fill_manual(values=colour)
    p <- p + xlab("dummy")
    p <- p + theme_bw()
    p <- p + theme(legend.title=element_blank())

    return(p)
}
    




                       
    
    
        

        


  
