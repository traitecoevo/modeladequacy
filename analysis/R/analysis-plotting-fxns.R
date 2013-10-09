## Generate plots for analysis
require(ggplot2)

cd <- getwd()
f <- file.path(cd, "output", "results-ml-sla-family.csv")
o <- file.path(cd, "output", "results-ml-sla-order.csv")
t <- file.path(cd, "output", "results-ml-sla-timeslice.csv")

fam <- read.csv(f, header=TRUE, as.is=TRUE)
ord <- read.csv(o, header=TRUE, as.is=TRUE)
tim <- read.csv(t, header=TRUE, as.is=TRUE)

res <- rbind(fam, ord, tim)



## little fxn to rename the rank for the time slice
clean.ml.results <- function(x){

    ## convert NAs to names so not discarded
    tmp <- which(is.na(x[,"rank"]))
    x[tmp, "rank"] <- rep("timeslice", length(tmp))

    tmp2 <- which(is.na(x[, "taxa"]))
    x[tmp2, "taxa"] <- rep("random", length(tmp2))

    ## remove all datapoints where mv.modelad
    x <- x[!is.na(x[,"mv.modelad"]), ]

    ## log mv.modelad
    x[,"mv.modelad"] <- log(x[,"mv.modelad"])

    ## make rank a factor
    x[,"rank"] <- as.factor(x[,"rank"])

    as.data.frame(x)
    
}



modelad.age.plot <- function(data){

    trait <- unique(data[,"trait"])
    
    .e <- environment()

    p <- ggplot(data, aes(log(age), mv.modelad), environment=.e)
    p <- p + geom_point(aes(colour=rank), size=3, alpha=0.9)
    p <- p + theme_bw()
    p <- p + xlab("Age of clade")
    p <- p + ylab("Log mahalanobis distance")
    p <- p + ggtitle(paste("Model adequacy versus clade age -", trait, sep=" "))
    print(p)
}




modelad.size.plot <- function(data){

    trait <- unique(data[,"trait"])
    
    .e <- environment()

    p <- ggplot(data, aes(log(size), mv.modelad), environment=.e)
    p <- p + geom_point(aes(colour=rank), size=3, alpha=0.9)
    p <- p + theme_bw()
    p <- p + xlab("Log number of taxa")
    p <- p + ylab("Log mahalanobis distance")
    p <- p + ggtitle(paste("Model adequacy versus clade size -", trait, sep=" "))
    print(p)
}



## Make plots
dd <- clean.ml.results(res)

## model adequacy versus age
modelad.age.plot(dd)

## model adequacy versus size
modelad.size.plot(dd)






    
