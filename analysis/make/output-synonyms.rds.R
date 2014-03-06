#!/usr/bin/env Rscript
library(methods)
dat <- read.csv("data/spermatophyta_synonyms_PLANTLIST.csv",
                stringsAsFactors=FALSE)
## Drop an extra column that has turned up for some reason:
dat <- dat[c("synonym", "species", "genus", "valid")]
saveRDS(dat, "output/synonyms.rds")
