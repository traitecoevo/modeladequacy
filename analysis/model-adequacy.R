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
lfn <- m$get("species_leaf_n")
sdm <- m$get("species_seed_mass")

fits_ml <- m$get("fits_ml")
fits_ml_best <- m$get("fits_ml_best")

fits_bayes <- m$get("fits_bayes")
fits_bayes_best <- m$get("fits_bayes_best")

e <- maker:::maker_environment(character(0), m)
maker_environment_attach(e)

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
## Overlap between tree and leaf_n
length(intersect(t$tip.label, lfn$gs))

## After splitting and fitting, number of clades in the dataset:
nrow(fits_ml)


## ## Results from case study: seed mass evolution in Meliaceae and Fagaceae

## Subset the ML analyses by dataset
## Meliaceae
mm <- fits_ml[fits_ml$taxa == "Meliaceae" & fits_ml$trait == "seed_mass", ]
## Fagaceae
ff <- fits_ml[fits_ml$taxa == "Fagaceae" & fits_ml$trait == "seed_mass", ]

## AIC weight of the OU model for Meliaceae
mm$aicw.ou

## AIC weight of the OU model for Fagaceae
ff$aicw.ou

## Get p-values for all six summary stats
pvalue_names <- paste0(pvalue_names_arbutus(), ".ml.ou")

## Meliaceae:
mm[, pvalue_names]

## Fagaceae:
ff[, pvalue_names]

## ## Results from conventional model comparison using AIC:

## Get AIC support for each model
aic <- build_ic(fits_ml)

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
aic_large <- aic[fits_ml$size >= 100, ]

## How many of the large trees supported by a single model with at
## least 90% AICw?
sum(apply(aic_large, 1, max) > 0.9)

## Of the large trees, how many support OU with at least 90% AICw
sum(aic_large$OU > 0.9)

## ### Code for plotting the relative AIC support for the different models

## Plot AIC model support
fig_model_support_ic(fits_ml)


## ## Results from conventional model comparison using DIC (Bayesian analysis):

## Get DIC support for each model
dic_names <- c("dicw.bm", "dicw.ou", "dicw.eb")
dic <- fits_bayes[, dic_names]
colnames(dic) <- c("BM", "OU", "EB")
dic <- dic[order(dic[,"OU"], dic[,"BM"], decreasing = TRUE),]

dic_best <- names(dic)[apply(dic, 1, which.max)]

## How many clades best supported by each model:
table(dic_best)

## How many clades best supported by OU
sum(dic_best == "OU")

## How many clades support OU with 100% of DICw
sum(dic$OU == 1)

## How many clades support OU with >75% of DICw
sum(dic$OU > 0.75)

## How many clades support EB with >75% of DICw
sum(dic$EB > 0.75)

## Now subsetting the taxa to only look at trees with at least 100 taxa
dic_large <- dic[fits_bayes$size >= 100,]

## How many of the large trees supported by a single model with at
## least 90% DICw
sum(apply(dic_large, 1, max) > 0.9)

## Of the large trees, how many support OU with at least 90% DICw
sum(dic_large$OU > 0.9)

## Relative support for different models using DIC
fig_model_support_dic(dic)


## ## Model adequacy results from fitting models with ML

## Use function prune_dataset_best_ml() to compile dataset of only the
## best supported of the three models for each dataset
ml_best <- prune_dataset_best(fits_ml)


## ### Overall model adequacy statistics

## For each dataset, how many summary statistics detected model violations
all_p_ml <- count_pvalues(ml_best, 0.05)

## How many datasets did not violate any of the summary statistics
sum(all_p_ml == 0)

## Number of datasets rejected by M_pic
sum(ml_best$m.pic <= 0.05)

## For those clades rejected by M_pic, in what proportion of cases was the
## value less than the median of the simulated values
mpic <- ml_best[ml_best$m.pic <= 0.05, ]
sum(mpic$mean.diag < 500)

## Number of datasets rejected by V_pic
sum(ml_best$v.pic <= 0.05)
## Number of datasets rejected by S.var
sum(ml_best$s.var <= 0.05)
## Number of datasets rejected by S.anc
sum(ml_best$s.anc <= 0.05)
## Number of datasets rejected by S.hgt
sum(ml_best$s.hgt <= 0.05)
## Number of datasets rejected by D.ks
sum(ml_best$d.ks <= 0.05)


## ### Model adequacy by trait

## Results from SLA data
ml_best_sla <- ml_best[ml_best$trait == "sla",]

## Number of datasets with SLA data
nrow(ml_best_sla)

## For each sla dataset, how many summary statistics detected model violations
## at p = 0.05
sla_p_ml <- count_pvalues(ml_best_sla, 0.05)
## How many rejected by at least one summary statistic
sum(sla_p_ml >= 1)
## How many rejected by at least two
sum(sla_p_ml >= 2)
## How many rejected by at least three
sum(sla_p_ml >= 3)

## Results from seed mass data
ml_best_sdm <- ml_best[ml_best$trait == "seed_size",]

## Number of datasets with seed mass data
nrow(ml_best_sdm)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
sdm_p_ml <- count_pvalues(ml_best_sdm, 0.05)
## How many rejected by at least one summary statistic
sum(sdm_p_ml >= 1)
## How many rejected by at least two
sum(sdm_p_ml >= 2)
## How many rejected by at least three
sum(sdm_p_ml >= 3)


## Results from leaf nitrogen data
ml_best_lfn <- ml_best[ml_best$trait == "leaf_n",]

## Number of datasets with leaf nitrogen data
nrow(ml_best_lfn)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
lfn_p_ml <- count_pvalues(ml_best_lfn, 0.05)
## How many rejected by at least one summary statistic
sum(lfn_p_ml >= 1)
## How many rejected by at least two
sum(lfn_p_ml >= 2)
## How many rejected by at least three
sum(lfn_p_ml >= 3)


## ### Code for plotting the distribution of p-values for all three traits

## Generate the plot
fig_pval_histogram(ml_best)



## ## Model adequacy results from fitting models with MCMC

## Use function prune.dataset.best.bayes() to compile a dataset of only the best
## supported model for each clade
bay.best <- prune.dataset.best.bayes(fits_bayes)


## ### Overall model adequacy statistics

## For each dataset, how many summary statistics detected model violations
all.p.bay <- count_pvalues(bay.best, 0.05)

## How many datasets did not violate any of the summary statistics
sum(all.p.bay == 0)

## Number of datasets rejected by M.pic
sum(bay.best$m.pic <= 0.05)


## Number of datasets rejected by V.pic
sum(bay.best$v.pic <= 0.05)
## Number of datasets rejected by S.var
sum(bay.best$s.var <= 0.05)
## Number of datasets rejected by S.anc
sum(bay.best$s.anc <= 0.05)
## Number of datasets rejected by S.hgt
sum(bay.best$s.hgt <= 0.05)
## Number of datasets rejected by D.ks
sum(bay.best$d.ks <= 0.05)


## ### Model adequacy by trait

## Results from SLA data
bay.best.sla <- bay.best[bay.best$trait == "sla",]

## For each sla dataset, how many summary statistics detected model violations
## at p = 0.05
sla.p.bay <- count_pvalues(bay.best.sla, 0.05)
## How many rejected by at least one summary statistic
sum(sla.p.bay >= 1)
## How many rejected by at least two
sum(sla.p.bay >= 2)
## How many rejected by at least three
sum(sla.p.bay >= 3)


## Results from seed mass data
bay.best.sdm <- bay.best[which(bay.best$trait == "seedmass"),]

## Number of datasets with seed mass data
nrow(bay.best.sdm)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
sdm.p.bay <- count_pvalues(bay.best.sdm, 0.05)
## How many rejected by at least one summary statistic
sum(sdm.p.bay >= 1)
## How many rejected by at least two
sum(sdm.p.bay >= 2)
## How many rejected by at least three
sum(sdm.p.bay >= 3)


## Results from leaf nitrogen data
bay.best.lfn <- bay.best[which(bay.best$trait == "leaf_n"),]

## Number of datasets with leaf nitrogen data
nrow(bay.best.lfn)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
lfn.p.bay <- count_pvalues(bay.best.lfn, 0.05)
## How many rejected by at least one summary statistic
sum(lfn.p.bay >= 1)
## How many rejected by at least two
sum(lfn.p.bay >= 2)
## How many rejected by at least three
sum(lfn.p.bay >= 3)


## Plotting the distribution of p-values for all three traits
fig_pval_histogram(bay.best)




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
aic_df <- build.table.adequacy.aic(fits_ml)

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
dic.df <- build.table.adequacy.dic(fits_bayes)

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
fig.modelad.size(ml_best)

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
fig.modelad.age(ml_best)

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
    
    to.pdf("output/figs/ad-aic.pdf", width=7, height=6,
           fig.modelad.aic(aic.df))

    to.pdf("output/figs/ad-dic.pdf", width=7, height=6,
           fig.modelad.dic(dic.df))

    to.pdf("output/figs/ad-size-ml.pdf", width=7, height=6,
           fig.modelad.size(ml_best))

    to.pdf("output/figs/ad-size-bayes.pdf", width=7, height=6,
           fig.modelad.size(bay.best))

    to.pdf("output/figs/ad-age-ml.pdf", width=7, height=6,
           fig.modelad.age(ml_best))

    to.pdf("output/figs/ad-age-bayes.pdf", width=7, height=6,
           fig.modelad.age(bay.best))

    if (file.info("output/results-ml")[,"isdir"])
        to.pdf("output/figs/two-clades.pdf", width=9, height=5,
               fig.two.clades())

}

