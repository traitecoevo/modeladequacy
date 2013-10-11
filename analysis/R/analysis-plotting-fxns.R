## Generate plots for analysis
require(ggplot2)

cd <- getwd()


## read in sla results
sla.fam <- read.csv(file.path(cd, "output", "results-ml-sla-family.csv"),
                header=TRUE, as.is=TRUE)
sla.ord <- read.csv(file.path(cd, "output", "results-ml-sla-order.csv"),
                header=TRUE, as.is=TRUE)
sla.time <- read.csv(file.path(cd, "output", "results-ml-sla-timeslice.csv"),
                header=TRUE, as.is=TRUE)

sla <- rbind(sla.fam, sla.ord, sla.time)

## read in seedMass results
sm.fam <- read.csv(file.path(cd, "output", "results-ml-seedMass-family.csv"),
                   header=TRUE, as.is=TRUE)
sm.ord <- read.csv(file.path(cd, "output", "results-ml-seedMass-order.csv"),
                   header=TRUE, as.is=TRUE)
sm.time <- read.csv(file.path(cd, "output", "results-ml-seedMass-timeslice.csv"),
                    header=TRUE, as.is=TRUE)

sm <- rbind(sm.fam, sm.ord, sm.time)





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

## Clean data
dat.sla <- clean.ml.results(sla)
dat.sm <- clean.ml.results(sm)

## model adequacy versus age
pdf(file.path(cd, "output", "results-ml-sla-adequacy-age.pdf"))
modelad.age.plot(dat.sla)
dev.off()


pdf(file.path(cd, "output", "results-ml-seedMass-adequacy-age.pdf"))
modelad.age.plot(dat.sm)
dev.off()


## model adequacy versus size
pdf(file.path(cd, "output", "results-ml-sla-adequacy-taxa.pdf"))
modelad.size.plot(dat.sla)
dev.off()


pdf(file.path(cd, "output", "results-ml-seedMass-adequacy-taxa.pdf"))
modelad.size.plot(dat.sm)
dev.off()







    
