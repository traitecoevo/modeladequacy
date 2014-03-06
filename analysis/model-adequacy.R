# Analysis of model adequacy for angiosperm functional traits
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

## Load in helper functions for analysis
source("R/model-adequacy-helper.R")

### Set options
##+ echo=FALSE, results=FALSE
knitr::opts_chunk$set(tidy=FALSE)

## Colours used throughout
col <- c("#a63813", "#4d697f", "gray15")


# Trait data across the angiosperm phylogeny
## import big tree
tree <- get.tree()

## extract angiosperms
t <- extract.clade(tree, node="Angiospermae")
## number of total angiosperms in tree
Ntip(t)

## age of tree
max(branching.times(t))

## read in the three data sets:
sla <- read.csv("output/species-mean-sla.csv", stringsAsFactors=FALSE)
sdm <- read.csv("output/species-mean-seedMass.csv", stringsAsFactors=FALSE)
lfn <- read.csv("output/species-mean-leafN.csv", stringsAsFactors=FALSE)

## number of species for which we have sla
nrow(sla)
## overlap between tree and sla
length(intersect(t$tip.label, sla$gs))

## number of species for which we have seedmass data
nrow(sdm)
## overlap between tree and seedmass
length(intersect(t$tip.label, sdm$gs))

## number of species for which we have leafn data
nrow(lfn)
## overlap between tree and leafN
length(intersect(t$tip.label, lfn$gs))


## TODO: Add function for making tree plots!!!




## Read in results from fitting models using maximum likelihood
ml <- read.csv("output/ml-results.csv", as.is=TRUE, row.names=1)

## Read in results from fitting models using MCMC
bay <- read.csv("output/bayes-results.csv", as.is=TRUE, row.names=1)

## Number of clades in dataset
nrow(ml)




# Results from conventional model comparison using AIC:

## Get AIC support for each model
aic.names <- c("aicw.bm", "aicw.ou", "aicw.eb")
aic <- ml[,aic.names]
colnames(aic) <- c("BM", "OU", "EB")
aic <- aic[order(aic[,"OU"], aic[,"BM"], decreasing = TRUE),]

aic.best <- sapply(seq_len(nrow(aic)), function(x)
                   return(colnames(aic)[which(aic[x,] == max(aic[x,]))]))

## How many clades best supported by OU
length(which(aic.best == "OU"))

## How many clades support OU with 100% of AICw
length(which(aic$OU == 1))

## How many clades support OU with >75% of AICw
length(which(aic$OU > 0.75))

## How many clades support EB with >75% of AICw
length(which(aic$EB > 0.75))

## Now subsetting the taxa to only look at trees with at least 100 taxa
lg.trees <- which(ml$size >= 100)
aic.lg <- aic[lg.trees,]

## How many of the large trees supported by a single model with at least 90% AICw
aic.lg.90 <- sapply(seq_len(nrow(aic.lg)), function(x)
                    return(max(aic.lg[x,]) > 0.9))
length(which(aic.lg.90))

## Of the large trees, how many support OU with at least 90% AICw
length(which(aic.lg$OU > 0.9))


# Plotting the relative AIC support for the different models
fig.model.support.aic <- function(aic){
    ## add dummy variable
    dd <- cbind(rownames(aic), aic)
    colnames(dd) <- c("Subclade", colnames(aic))
    ord <- dd[order(dd[,"OU"], dd[,"BM"], decreasing = TRUE),"Subclade"]
    ## reorient data frame for plotting
    df <- melt(dd)
    colnames(df)[2:3] <- c("model", "weight")
    df$Subclade <- factor(df$Subclade, ord)

    ## create geom_bar plot
    .e <- environment()
    p <- ggplot(df, aes(factor(Subclade), weight, fill=model, order=model), environment=.e)
    p <- p + geom_bar(stat="identity", position="stack",width=1)
    p <- p + scale_y_continuous(name="AIC weight")
    p <- p + scale_fill_manual(values=col)
    p <- p + theme_bw()
    p <- p + theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.y=element_blank(),
                   strip.background=element_rect(fill="white"),
                   plot.background=element_blank())
    p <- p + xlab("Dataset")
    p
}

fig.model.support.aic(aic)





# Results from conventional model comparison using DIC (Bayesian analysis):

## Get DIC support for each model
dic.names <- c("dicw.bm", "dicw.ou", "dicw.eb")
dic <- bay[,dic.names]
colnames(dic) <- c("BM", "OU", "EB")
dic <- dic[order(dic[,"OU"], dic[,"BM"], decreasing = TRUE),]

dic.best <- sapply(seq_len(nrow(dic)), function(x)
                   return(colnames(dic)[which(dic[x,] == max(dic[x,]))]))

## How many clades best supported by OU
length(which(dic.best == "OU"))

## How many clades support OU with 100% of DICw
length(which(dic$OU == 1))

## How many clades support OU with >75% of AICw
length(which(dic$OU > 0.75))

## How many clades support EB with >75% of AICw
length(which(dic$EB > 0.75))

## Now subsetting the taxa to only look at trees with at least 100 taxa
lg.trees.bay <- which(bay$size >= 100)
dic.lg <- dic[lg.trees.bay,]

## How many of the large trees supported by a single model with at least 90% DICw
dic.lg.90 <- sapply(seq_len(nrow(dic.lg)), function(x)
                    return(max(dic.lg[x,]) > 0.9))
length(which(dic.lg.90))

## Of the large trees, how many support OU with at least 90% DICw
length(which(dic.lg$OU > 0.9))


# Plotting the relative DIC support for the different models
fig.model.support.dic <- function(dic){
    ## add dummy variable
    dd <- cbind(rownames(dic), dic)
    colnames(dd) <- c("Subclade", colnames(dic))
    ord <- dd[order(dd[,"OU"], dd[,"BM"], decreasing = TRUE),"Subclade"]
    ## reorient data frame for plotting
    df <- melt(dd)
    colnames(df)[2:3] <- c("model", "weight")
    df$Subclade <- factor(df$Subclade, ord)

    ## create geom_bar plot
    .e <- environment()
    p <- ggplot(df, aes(factor(Subclade), weight, fill=model, order=model), environment=.e)
    p <- p + geom_bar(stat="identity", position="stack",width=1)
    p <- p + scale_y_continuous(name="DIC weight")
    p <- p + scale_fill_manual(values=col)
    p <- p + theme_bw()
    p <- p + theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.y=element_blank(),
                   strip.background=element_rect(fill="white"),
                   plot.background=element_blank())
    p <- p + xlab("Dataset")
    p
}

fig.model.support.dic(dic)







# Model adequacy results from fitting models with ML

## Only use the best supported model
ml.best <- prune.dataset.best.ml(ml)


## Across all three traits

## For each dataset, how many summary statistics detected model violations
all.p.ml <- count.pvalues(ml.best, 0.05)

## How many datasets did not violate any of the summary statistics
length(which(all.p.ml == 0))

## Number of datasets rejected by M_pic
length(which(ml.best$m.pic <= 0.05))

## For those clades rejected by M_pic, in what proportion of cases was the
## value less than the median of the simulated values
mpic <- ml.best[which(ml.best$m.pic <= 0.05), ]
length(which(mpic$mean.diag < 500))

## Number of datasets rejected by V_pic
length(which(ml.best$v.pic <= 0.05))
## Number of datasets rejected by S.var
length(which(ml.best$s.var <= 0.05))
## Number of datasets rejected by S.anc
length(which(ml.best$s.anc <= 0.05))
## Number of datasets rejected by S.hgt
length(which(ml.best$s.hgt <= 0.05))
## Number of datasets rejected by D.ks
length(which(ml.best$d.ks <= 0.05))


## Looking at each trait independently

## Results from SLA data
ml.best.sla <- ml.best[which(ml.best$trait == "SLA"),]

## Number of datasets with SLA data
nrow(ml.best.sla)

## For each sla dataset, how many summary statistics detected model violations
## at p = 0.05
sla.p.ml <- count.pvalues(ml.best.sla, 0.05)
## How many rejected by at least one summary statistic
length(which(sla.p.ml >= 1))
## How many rejected by at least two
length(which(sla.p.ml >= 2))
## How many rejected by at least three
length(which(sla.p.ml >= 3))


## Results from seed mass data
ml.best.sdm <- ml.best[which(ml.best$trait == "seedmass"),]

## Number of datasets with seed mass data
nrow(ml.best.sdm)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
sdm.p.ml <- count.pvalues(ml.best.sdm, 0.05)
## How many rejected by at least one summary statistic
length(which(sdm.p.ml >= 1))
## How many rejected by at least two
length(which(sdm.p.ml >= 2))
## How many rejected by at least three
length(which(sdm.p.ml >= 3))


## Results from leaf nitrogen data
ml.best.lfn <- ml.best[which(ml.best$trait == "leafn"),]

## Number of datasets with leaf nitrogen data
nrow(ml.best.lfn)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
lfn.p.ml <- count.pvalues(ml.best.lfn, 0.05)
## How many rejected by at least one summary statistic
length(which(lfn.p.ml >= 1))
## How many rejected by at least two
length(which(lfn.p.ml >= 2))
## How many rejected by at least three
length(which(lfn.p.ml >= 3))


# Plotting the distribution of p-values for all three traits
fig.pval.histogram <- function(best){
    ## prune out irrelevant categories
    best <- cbind(rownames(best), best[,c("trait", "m.pic", "v.pic", "s.var",
                                    "s.anc", "s.hgt", "d.ks")])

 
    ## melt
    df <- melt(best)

    ## set the order and the labels
    df$trait <- factor(df$trait, levels=c("SLA", "seedmass", "leafn"),
                       labels=c("SLA", "SeedMass", "LeafN"))
                        
    df$variable <- factor(df$variable,
                          labels=c("italic(M[PIC])", "italic(V[PIC])", "italic(S[VAR])", "italic(S[ANC])", "italic(S[HGT])", "italic(D[KS])"))
    .e <- environment()

    p <- ggplot(df, aes(x=value), environment = .e)
    p <- p + geom_histogram(binwidth=0.025, alpha=0.8, aes(y=..density.., fill=factor(trait)))
    p <- p + scale_fill_manual(values=col)
    p <- p + theme_bw()
    p <- p + xlab("p-value")
    p <- p + facet_grid(trait~variable, labeller = label_parsed)
    p <- p + theme(strip.background=element_rect(fill="white"),
                   plot.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   axis.text=element_text(size=8),
                   axis.ticks.y=element_blank(),
                   axis.text.y=element_blank(),
                   legend.position="none")
    p
}

fig.pval.histogram(ml.best)









# Model adequacy results from fitting models with MCMC

## Only use the best supported model
bay.best <- prune.dataset.best.bayes(bay)


## Across all three traits

## For each dataset, how many summary statistics detected model violations
all.p.bay <- count.pvalues(bay.best, 0.05)

## How many datasets did not violate any of the summary statistics
length(which(all.p.bay == 0))

## Number of datasets rejected by M_pic
length(which(bay.best$m.pic <= 0.05))


## Number of datasets rejected by V_pic
length(which(bay.best$v.pic <= 0.05))
## Number of datasets rejected by S.var
length(which(bay.best$s.var <= 0.05))
## Number of datasets rejected by S.anc
length(which(bay.best$s.anc <= 0.05))
## Number of datasets rejected by S.hgt
length(which(bay.best$s.hgt <= 0.05))
## Number of datasets rejected by D.ks
length(which(bay.best$d.ks <= 0.05))


## Looking at each trait independently

## Results from SLA data
bay.best.sla <- bay.best[which(bay.best$trait == "SLA"),]

## For each sla dataset, how many summary statistics detected model violations
## at p = 0.05
sla.p.bay <- count.pvalues(bay.best.sla, 0.05)
## How many rejected by at least one summary statistic
length(which(sla.p.bay >= 1))
## How many rejected by at least two
length(which(sla.p.bay >= 2))
## How many rejected by at least three
length(which(sla.p.bay >= 3))


## Results from seed mass data
bay.best.sdm <- ml.best[which(bay.best$trait == "seedmass"),]

## Number of datasets with seed mass data
nrow(bay.best.sdm)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
sdm.p.bay <- count.pvalues(bay.best.sdm, 0.05)
## How many rejected by at least one summary statistic
length(which(sdm.p.bay >= 1))
## How many rejected by at least two
length(which(sdm.p.bay >= 2))
## How many rejected by at least three
length(which(sdm.p.bay >= 3))


## Results from leaf nitrogen data
bay.best.lfn <- bay.best[which(bay.best$trait == "leafn"),]

## Number of datasets with leaf nitrogen data
nrow(bay.best.lfn)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
lfn.p.bay <- count.pvalues(bay.best.lfn, 0.05)
## How many rejected by at least one summary statistic
length(which(lfn.p.bay >= 1))
## How many rejected by at least two
length(which(lfn.p.bay >= 2))
## How many rejected by at least three
length(which(lfn.p.bay >= 3))


# Plotting the distribution of p-values for all three traits
fig.pval.histogram(bay.best)





# Model adequacy versus size
## Plot a multivariate measure of model adequacy (Mahalanobis distance) against clade size

fig.modelad.size <- function(df){

    ## Capitalize ranks
    df$rank <- sapply(df$rank, function(x) cap.ranks(x))
    ## rename and reorder trait
    df$trait <- sapply(df$trait, function(x) rename.traits(x))
    df$trait <- factor(df$trait, c("SLA", "SeedMass", "LeafN"))
    
   .e <- environment()

    ## the occasional dataset may have NA for Mahalanobis distance
    ## remove this for the plot
    df <- na.omit(df)

    p <- ggplot(df, aes(size, mv), environment=.e)
    p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.6)

    p <- p + scale_colour_manual(values=col)
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   axis.line = element_line(color = 'black'))
    p <- p + scale_y_log10()
    p <- p + scale_x_log10()
    p <- p + xlab("Number of taxa")
    p <- p + ylab("Mahalanobis distance")
    p
}


## For the best supported model from ML analysis
fig.modelad.size(ml.best)

## For the best supported model from the Bayesian analysis
fig.modelad.size(bay.best)


# Model adequacy versus age
## Plot a multivariate measure of model adequacy (Mahalanobis distance) against clade age

fig.modelad.age <- function(df){

    ## Capitalize ranks
    df$rank <- sapply(df$rank, function(x) cap.ranks(x))
    ## rename and reorder trait
    df$trait <- sapply(df$trait, function(x) rename.traits(x))
    df$trait <- factor(df$trait, c("SLA", "SeedMass", "LeafN"))
    
   .e <- environment()

    ## the occasional dataset may have NA for Mahalanobis distance
    ## remove this for the plot
    df <- na.omit(df)

    p <- ggplot(df, aes(age, mv), environment=.e)
    p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.6)

    p <- p + scale_colour_manual(values=col)
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   axis.line = element_line(color = 'black'))
    p <- p + scale_y_log10()
    p <- p + scale_x_log10()
    p <- p + xlab("Age of crown group (my)")
    p <- p + ylab("Mahalanobis distance")
    p
}


## For the best supported model from ML analysis
fig.modelad.age(ml.best)

## For the best supported model from the Bayesian analysis
fig.modelad.age(bay.best)


## produce figures
## TODO: THIS CURRENTLY DOES WORK TO PRODUCE PDFs. I AM PROBABY MISSING A FUNCTION
if (!interactive()){

    ## create a directory for the figures
    dir.create("output/figs", FALSE)
    
    to.pdf("output/figs/aic-support.pdf", width=9, height=6,
           fig.model.support.aic(aic))

    to.pdf("output/figs/dic-support.pdf", width=9, height=6,
           fig.model.support.dic(dic))

    to.pdf("output/figs/pval-hist-ml.pdf", width=9, height=5,
           fig.pval.histogram(ml.best), onefile=FALSE)

    to.pdf("output/figs/pval-hist-bayes.pdf", width=9, height=5,
           fig.pval.histogram(bay.best), onefile=FALSE)

    to.pdf("output/figs/ad-size-ml.pdf", width=7, height=6,
           fig.modelad.size(ml.best))

    to.pdf("output/figs/ad-size-bayes.pdf", width=7, height=6,
           fig.modelad.size(bay.best))

    to.pdf("output/figs/ad-age-ml.pdf", width=7, height=6,
           fig.modelad.age(ml.best))

    to.pdf("output/figs/ad-age-bayes.pdf", width=7, height=6,
           fig.modelad.age(bay.best))

}

