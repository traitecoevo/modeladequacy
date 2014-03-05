#!/usr/bin/env Rscript
source("R/read-data-functions.R")
saveRDS(build.data("seedMass"), "output/data-seedMass.rds")
