## plotting functions for arbutus summary statistics



## plot slope of contrasts versus variance
## not really that flexible at the moment

contrastVar.plot <- function(x, col=c("dodgerblue4", "darkblue"), ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## create data frame for pics
    dat <- as.data.frame(ut$pics)
    dat$contrasts <- abs(dat$contrasts)

    p <- ggplot(dat, aes(variance, contrasts))
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

    p <- ggplot(dat, aes(asr, pics))
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

    p <- ggplot(dat, aes(nh, pics))
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

