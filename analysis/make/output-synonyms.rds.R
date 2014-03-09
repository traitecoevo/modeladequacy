#!/usr/bin/env Rscript
library(methods)
source("R/read-data-functions.R")
tree <- get.tree()

dat <- read.csv("data/spermatophyta_synonyms_PLANTLIST.csv",
                stringsAsFactors=FALSE)
## Drop an extra column that has turned up for some reason:
dat <- dat[c("synonym", "species", "genus")]
dat <- dat[dat$synonym != dat$species,]
dat$synonym[dat$synonym %in% tree$tip.label] <- NA
saveRDS(dat, "output/synonyms.rds")
