## # Model adequacy and the macroevolution of angiosperm functional traits -- Analysis

## Load packages and helper functions
## source("R/model-adequacy-analysis.R")

### Set options
##+ echo=FALSE, results=FALSE
knitr::opts_chunk$set(tidy=FALSE)

## ## Trait data across the angiosperm phylogeny

## Extract angiosperms from the phylogeny
t <- extract.clade(vascular_plant_phylogeny, node="Angiospermae")

## Number of taxa in angiosperm tree
Ntip(t)

## Age of tree
max(branching.times(t))

## We have three species-level data sets: specific leaf area (metres
## square per gram; sla), seed mass (gram; sdm), and leaf nitrogen (grams per
## gram; lfn).

sla <- species_sla
lfn <- species_leaf_n
sdm <- species_seed_mass

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
dic <- build_ic(fits_bayes)

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
fig_model_support_ic(fits_bayes)

## ## Model adequacy results from fitting models with ML

## Use function prune_dataset_best_ml() to compile dataset of only the
## best supported of the three models for each dataset

## ### Overall model adequacy statistics

## For each dataset, how many summary statistics detected model violations
all_p_ml <- count_pvalues(fits_ml_best, 0.05)

## How many datasets did not violate any of the summary statistics
sum(all_p_ml == 0)

## Number of datasets rejected by M.sig
sum(fits_ml_best$m.sig <= 0.05)

## For those clades rejected by M.sig, in what proportion of cases was the
## value less than the median of the simulated values
msig <- fits_ml_best[fits_ml_best$m.sig <= 0.05, ]
sum(msig$mean.diag < 500)

## Number of datasets rejected by C.var
sum(fits_ml_best$c.var <= 0.05)
## Number of datasets rejected by S.var
sum(fits_ml_best$s.var <= 0.05)
## Number of datasets rejected by S.asr
sum(fits_ml_best$s.asr <= 0.05)
## Number of datasets rejected by S.hgt
sum(fits_ml_best$s.hgt <= 0.05)
## Number of datasets rejected by D.cdf
sum(fits_ml_best$d.cdf <= 0.05)


## ### Model adequacy by trait

## Results from SLA data
fits_ml_best_sla <- fits_ml_best[fits_ml_best$trait == "sla",]

## Number of datasets with SLA data
nrow(fits_ml_best_sla)

## For each sla dataset, how many summary statistics detected model violations
## at p = 0.05
sla_p_ml <- count_pvalues(fits_ml_best_sla, 0.05)
## How many rejected by at least one summary statistic
sum(sla_p_ml >= 1)
## How many rejected by at least two
sum(sla_p_ml >= 2)
## How many rejected by at least three
sum(sla_p_ml >= 3)

## Results from seed mass data
fits_ml_best_sdm <- fits_ml_best[fits_ml_best$trait == "seed_mass",]

## Number of datasets with seed mass data
nrow(fits_ml_best_sdm)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
sdm_p_ml <- count_pvalues(fits_ml_best_sdm, 0.05)
## How many rejected by at least one summary statistic
sum(sdm_p_ml >= 1)
## How many rejected by at least two
sum(sdm_p_ml >= 2)
## How many rejected by at least three
sum(sdm_p_ml >= 3)


## Results from leaf nitrogen data
fits_ml_best_lfn <- fits_ml_best[fits_ml_best$trait == "leaf_n",]

## Number of datasets with leaf nitrogen data
nrow(fits_ml_best_lfn)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
lfn_p_ml <- count_pvalues(fits_ml_best_lfn, 0.05)
## How many rejected by at least one summary statistic
sum(lfn_p_ml >= 1)
## How many rejected by at least two
sum(lfn_p_ml >= 2)
## How many rejected by at least three
sum(lfn_p_ml >= 3)


## ### Code for plotting the distribution of p-values for all three traits

## Generate the plot
fig_pval_histogram(fits_ml_best)



## ## Model adequacy results from fitting models with MCMC

## Use function prune.dataset.best.bayes() to compile a dataset of only the best
## supported model for each clade


## ### Overall model adequacy statistics

## For each dataset, how many summary statistics detected model violations
all_p_bay <- count_pvalues(fits_bayes_best, 0.05)

## How many datasets did not violate any of the summary statistics
sum(all_p_bay == 0)

## Number of datasets rejected by M.sig
sum(fits_bayes_best$m.sig <= 0.05)


## Number of datasets rejected by C.var
sum(fits_bayes_best$c.var <= 0.05)
## Number of datasets rejected by S.var
sum(fits_bayes_best$s.var <= 0.05)
## Number of datasets rejected by S.asr
sum(fits_bayes_best$s.asr <= 0.05)
## Number of datasets rejected by S.hgt
sum(fits_bayes_best$s.hgt <= 0.05)
## Number of datasets rejected by D.cdf
sum(fits_bayes_best$d.cdf <= 0.05)


## ### Model adequacy by trait

## Results from SLA data
fits_bayes_best_sla <- fits_bayes_best[fits_bayes_best$trait == "sla",]

## For each sla dataset, how many summary statistics detected model violations
## at p = 0.05
sla_p_bay <- count_pvalues(fits_bayes_best_sla, 0.05)
## How many rejected by at least one summary statistic
sum(sla_p_bay >= 1)
## How many rejected by at least two
sum(sla_p_bay >= 2)
## How many rejected by at least three
sum(sla_p_bay >= 3)


## Results from seed mass data
fits_bayes_best_sdm <- fits_bayes_best[fits_bayes_best$trait == "seed_mass",]

## Number of datasets with seed mass data
nrow(fits_bayes_best_sdm)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
sdm_p_bay <- count_pvalues(fits_bayes_best_sdm, 0.05)
## How many rejected by at least one summary statistic
sum(sdm_p_bay >= 1)
## How many rejected by at least two
sum(sdm_p_bay >= 2)
## How many rejected by at least three
sum(sdm_p_bay >= 3)


## Results from leaf nitrogen data
fits_bayes_best_lfn <- fits_bayes_best[which(fits_bayes_best$trait == "leaf_n"),]

## Number of datasets with leaf nitrogen data
nrow(fits_bayes_best_lfn)

## For each seed mass dataset, how many summary statistics detected model violations
## at p = 0.05
lfn_p_bay <- count_pvalues(fits_bayes_best_lfn, 0.05)
## How many rejected by at least one summary statistic
sum(lfn_p_bay >= 1)
## How many rejected by at least two
sum(lfn_p_bay >= 2)
## How many rejected by at least three
sum(lfn_p_bay >= 3)


## Plotting the distribution of p-values for all three traits
fig_pval_histogram(fits_bayes_best)

## Compile table with the best model only, excluding datasets where BM
## is the best supported model. Calculate difference in AIC score
## between best supported alternate model and BM, plotting this
## against Mahalanobis distance:
fig_modelad_ic(fits_ml)

## As above, for DIC:
fig_modelad_ic(fits_bayes)

## ## Model adequacy versus size

## For the best supported model from ML analysis
fig_modelad_size(fits_ml_best)

## For the best supported model from the Bayesian analysis
fig_modelad_size(fits_bayes_best)

## ## Model adequacy versus age

## For the best supported model from ML analysis
fig_modelad_age(fits_ml_best)

## For the best supported model from the Bayesian analysis
fig_modelad_age(fits_bayes_best)

## ## Two clade example figure
fig_two_clades(example_fits)

