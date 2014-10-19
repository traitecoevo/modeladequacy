## # Model adequacy and the macroevolution of angiosperm functional traits -- Analysis

## Load packages and helper functions
source("R/model-adequacy-analysis.R")

## This is temporary:
source("maker/figures.R")

### Set options
##+ echo=FALSE, results=FALSE
knitr::opts_chunk$set(tidy=FALSE)

## Define colours used throughout
col <- c("#F46D43", "#3288BD", "#CDCD00")

## Haul in stuff from maker, at least for now:
library(maker)
m <- maker$new()
tree <- m$get("vascular_plant_phylogeny")

sla <- m$get("species_sla")
lfn <- m$get("species_leafN")
sdm <- m$get("species_seed_mass")

ml <- m$get("fits_ml")

## ## Trait data across the angiosperm phylogeny

## Extract angiosperms from the phylogeny
t <- extract.clade(tree, node="Angiospermae")

## Number of taxa in angiosperm tree
Ntip(t)

## Age of tree
max(branching.times(t))

## We have three species-level data sets: specific leaf area (metres
## square per gram; sla), seed mass (gram; sdm), and leaf nitrogen (grams per
## gram; lfn).

## Number of species for which we have sla
nrow(sla)
## Overlap between tree and sla
length(intersect(t$tip.label, sla$gs))

## Number of species for which we have seedmass data
nrow(sdm)
## Overlap between tree and seedmass
length(intersect(t$tip.label, sdm$gs))

## Number of species for which we have leafn data
nrow(lfn)
## Overlap between tree and leafN
length(intersect(t$tip.label, lfn$gs))

## After splitting and fitting, number of clades in the dataset:
nrow(ml)


## ## Results from case study: seed mass evolution in Meliaceae and Fagaceae

## Subset the ML analyses by dataset
## Meliaceae
mm <- ml[ml$taxa == "Meliaceae" & ml$trait == "seed_size", ]
## Fagaceae
ff <- ml[ml$taxa == "Fagaceae" & ml$trait == "seed_size", ]

## AIC weight of the OU model for Meliaceae
mm$aicw.ou

## AIC weight of the OU model for Fagaceae
ff$aicw.ou

## Get p-values for all six summary stats
pvalue_names_arbutus <- c("m.sig", "c.var", "s.var", "s.asr", "s.hgt", "d.cdf")
pvalue_names <- paste0(pvalue_names_arbutus, ".ml.ou")

## Meliaceae:
mm[, pvalue_names]

## Fagaceae:
ff[, pvalue_names]

## ## Results from conventional model comparison using AIC:

## Get AIC support for each model
aic_names <- c("aicw.bm", "aicw.ou", "aicw.eb")
aic <- ml[aic_names]
colnames(aic) <- c("BM", "OU", "EB")
aic <- aic[order(aic[,"OU"], aic[,"BM"], decreasing = TRUE),]

aic_best <- names(aic)[apply(aic, 1, which.max)]

## How many clades best supported by each model:
table(aic_best)

## How many clades support OU with 100% of AICw
sum(aic$OU == 1)

## How many clades support OU with >75% of AICw
sum(aic$OU > 0.75)

## How many clades support EB with >75% of AICw
sum(aic$EB > 0.75)

## Now subsetting the taxa to only look at trees with at least 100 taxa
aic_large <- aic[ml$size >= 100, ]

## How many of the large trees supported by a single model with at
## least 90% AICw?
sum(apply(aic_large, 1, max) > 0.9)

## Of the large trees, how many support OU with at least 90% AICw
sum(aic_large$OU > 0.9)

## ### Code for plotting the relative AIC support for the different models

## Plot AIC model support
fig_model_support_aic(aic)





## ## Results from conventional model comparison using DIC (Bayesian analysis):

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

## How many clades support OU with >75% of DICw
length(which(dic$OU > 0.75))

## How many clades support EB with >75% of DICw
length(which(dic$EB > 0.75))

## Now subsetting the taxa to only look at trees with at least 100 taxa
dic.lg <- dic[bay$size >= 100,]

## How many of the large trees supported by a single model with at least 90% DICw
dic.lg.90 <- sapply(seq_len(nrow(dic.lg)), function(x)
                    return(max(dic.lg[x,]) > 0.9))
length(which(dic.lg.90))

## Of the large trees, how many support OU with at least 90% DICw
length(which(dic.lg$OU > 0.9))


## ### Code for plotting the relative DIC support for the different models
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
    p <- p + scale_fill_manual("Model", values=col)
    p <- p + theme_bw()
    p <- p + theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.y=element_blank(),
                   panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.border=element_blank(),
                   panel.background=element_blank(),
                   strip.background=element_rect(fill="white"),
                   plot.background=element_blank())
    p <- p + xlab("Dataset")
    p
}

## Plot DIC model support
fig.model.support.dic(dic)



## ## Model adequacy results from fitting models with ML

## Use function prune.dataset.best.ml() to compile dataset of only the best supported of
## the three models for each dataset
ml.best <- prune.dataset.best.ml(ml)


## ### Overall model adequacy statistics

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


## ### Model adequacy by trait

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


## ### Code for plotting the distribution of p-values for all three traits
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
                          labels=c("italic(M[SIG])", "italic(C[VAR])", "italic(S[VAR])", "italic(S[ASR])", "italic(S[HGT])", "italic(D[CDF])"))
    .e <- environment()

    p <- ggplot(df, aes(x=value), environment = .e)
    p <- p + geom_histogram(binwidth=0.025, alpha=0.8, aes(y=..density.., fill=factor(trait)))
    p <- p + scale_fill_manual(values=col)
    p <- p + theme_bw()
    p <- p + xlab("p-value")
    p <- p + ylab("Density")
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

## Generate the plot
fig.pval.histogram(ml.best)



## ## Model adequacy results from fitting models with MCMC

## Use function prune.dataset.best.bayes() to compile a dataset of only the best
## supported model for each clade
bay.best <- prune.dataset.best.bayes(bay)


## ### Overall model adequacy statistics

## For each dataset, how many summary statistics detected model violations
all.p.bay <- count.pvalues(bay.best, 0.05)

## How many datasets did not violate any of the summary statistics
length(which(all.p.bay == 0))

## Number of datasets rejected by M.pic
length(which(bay.best$m.pic <= 0.05))


## Number of datasets rejected by V.pic
length(which(bay.best$v.pic <= 0.05))
## Number of datasets rejected by S.var
length(which(bay.best$s.var <= 0.05))
## Number of datasets rejected by S.anc
length(which(bay.best$s.anc <= 0.05))
## Number of datasets rejected by S.hgt
length(which(bay.best$s.hgt <= 0.05))
## Number of datasets rejected by D.ks
length(which(bay.best$d.ks <= 0.05))


## ### Model adequacy by trait

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
bay.best.sdm <- bay.best[which(bay.best$trait == "seedmass"),]

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


## Plotting the distribution of p-values for all three traits
fig.pval.histogram(bay.best)




## ### Code for plotting the relative support (AIC) for the best model (compared to BM) vs. a multivariate measure of model adequacy
fig.modelad.aic <- function(df){
    df <- prepare.df.for.ggplot(df)
    
   .e <- environment()

    ## need to set options for scientific notation
    options(scipen=1000)
    
    ## the occasional dataset may have NA for Mahalanobis distance
    ## remove this for the plot
    df <- na.omit(df)

    p <- ggplot(df, aes(diff.bm, mv), environment=.e)
    p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.8)
    p <- p + scale_colour_manual("Trait", values=col)
    p <- p + scale_shape_manual("Rank", values=c(15,16,17))
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   legend.justification=c(0.02,0.98),
                   legend.position=c(0.02,0.98),
                   legend.key=element_blank(),
                   axis.line = element_line(color = 'black'))
    p <- p + scale_y_log10()
    p <- p + scale_x_log10()
    p <- p + xlab("AIC(BM) - AIC(OU/EB)")
    p <- p + ylab("Mahalanobis distance")
    p
}

## Compile table with the best model only, excluding datasets where BM is the
## best supported model. Calculate difference in AIC score between best supported
## alternate model and BM.
aic_df <- build.table.adequacy.aic(ml)

fig.modelad.aic(aic_df)



## ### Code for plotting the relative support (DIC) for the best model (compared to BM) vs. a multivariate measure of model adequacy
fig.modelad.dic <- function(df){
    df <- prepare.df.for.ggplot(df)
    .e <- environment()

    ## the occasional dataset may have NA for Mahalanobis distance
    ## remove this for the plot
    df <- na.omit(df)

    p <- ggplot(df, aes(diff.bm, mv), environment=.e)
    p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.8)
    p <- p + scale_colour_manual("Trait", values=col)
    p <- p + scale_shape_manual("Rank", values=c(15,16,17))
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   legend.justification=c(0.02,0.98),
                   legend.position=c(0.02,0.98),
                   legend.key=element_blank(),
                   axis.line = element_line(color = 'black'))
    p <- p + scale_y_log10()
    p <- p + scale_x_log10()
    p <- p + xlab("DIC(BM) - DIC(OU/EB)")
    p <- p + ylab("Mahalanobis distance")
    p
}

## Compile table with the best model only, excluding datasets where BM is the
## best supported model. Calculate difference in DIC score between best supported
## alternate model and BM.
dic.df <- build.table.adequacy.dic(bay)

fig.modelad.dic(dic.df)







## ## Model adequacy versus size

## ### Code to plot a multivariate measure of model adequacy (Mahalanobis distance) against clade size
fig.modelad.size <- function(df){
    df <- prepare.df.for.ggplot(df)

   .e <- environment()

    ## the occasional dataset may have NA for Mahalanobis distance
    ## remove this for the plot
    df <- na.omit(df)

    p <- ggplot(df, aes(size, mv), environment=.e)
    p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.8)

    p <- p + scale_colour_manual("Trait", values=col)
    p <- p + scale_shape_manual("Rank", values=c(15,16,17))
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   legend.justification=c(0.98,0.02),
                   legend.position=c(0.98,0.02),
                   legend.key=element_blank(),
                   axis.line = element_line(color = 'black'))
    p <- p + scale_y_log10()
    p <- p + scale_x_log10()
    p <- p + xlab("Number of taxa")
    p <- p + ylab("Mahalanobis distance between observed and simulated")
    p
}


## For the best supported model from ML analysis
fig.modelad.size(ml.best)

## For the best supported model from the Bayesian analysis
fig.modelad.size(bay.best)


## ## Model adequacy versus age

## ### Code to plot a multivariate measure of model adequacy (Mahalanobis distance) against clade age
fig.modelad.age <- function(df){
    df <- prepare.df.for.ggplot(df)
    
   .e <- environment()

    ## the occasional dataset may have NA for Mahalanobis distance
    ## remove this for the plot
    df <- na.omit(df)

    p <- ggplot(df, aes(age, mv), environment=.e)
    p <- p + geom_point(aes(colour=trait, shape=rank), size=3, alpha=0.8)
    p <- p + scale_colour_manual("Trait", values=col)
    p <- p + scale_shape_manual("Rank", values=c(15,16,17))
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   legend.justification=c(0.02,0.98),
                   legend.position=c(0.02,0.98),
                   legend.key=element_blank(),
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


## ## Two clade example figure

## Only do this if results have actually been run
if (file.info("output/results-ml")[,"isdir"]){

## Read in Meliaceae and Fagaceae data
fig.two.clades <- function(){
    me.dat <- readRDS("output/results-ml/seedmass_clade_Meliaceae.rds")$OU
    fa.dat <- readRDS("output/results-ml/seedmass_clade_Fagaceae.rds")$OU
    par(mfrow=c(2,6))
    lapply(colnames(me.dat$ma$obs), function(x){
    par(mar=c(4,1,1,1))

    profiles.plot(me.dat$ma$sim[x], col.line=col[1],
                  opacity = 0.9, frame.plot=FALSE, yaxt="n",
                  xlab="", ylab="");
    abline(v=me.dat$ma$obs[,x], lty=2, lwd=2, col=col[2])})

    lapply(colnames(fa.dat$ma$obs), function(x){
    par(mar=c(4.5,1,1,1))
    xl <- strsplit(toupper(x), split=".")
    xlab <- bquote(.(xl[[1]][1], xl[[1]][2]) ~ [x
    
    if (x == "m.sig"){
        profiles.plot(fa.dat$ma$sim[x], col.line=col[3],
                  opacity = 0.9, frame.plot=FALSE, yaxt="n",
                  xlab=expression(italic(xl.f)[xl.s]), ylab="", cex.lab=1.5,
                  xlim=c(as.numeric(fa.dat$ma$obs[x]-0.05), max(fa.dat$ma$sim[x])))
    } else {
        profiles.plot(fa.dat$ma$sim[x], col.line=col[3],
                  opacity = 0.9, frame.plot=FALSE, yaxt="n",
                  xlab=expression(italic(xl.f)[xl.s]), ylab="", cex.lab=1.5)
    }
    abline(v=fa.dat$ma$obs[,x], lty=2, lwd=2, col=col[2])
    })
}

}





## ## Produce figures

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

    to.pdf("output/figs/ad-aic.pdf", width=7, height=6,
           fig.modelad.aic(aic.df))

    to.pdf("output/figs/ad-dic.pdf", width=7, height=6,
           fig.modelad.dic(dic.df))

    to.pdf("output/figs/ad-size-ml.pdf", width=7, height=6,
           fig.modelad.size(ml.best))

    to.pdf("output/figs/ad-size-bayes.pdf", width=7, height=6,
           fig.modelad.size(bay.best))

    to.pdf("output/figs/ad-age-ml.pdf", width=7, height=6,
           fig.modelad.age(ml.best))

    to.pdf("output/figs/ad-age-bayes.pdf", width=7, height=6,
           fig.modelad.age(bay.best))

    if (file.info("output/results-ml")[,"isdir"])
        to.pdf("output/figs/two-clades.pdf", width=9, height=5,
               fig.two.clades())

}

