#!/usr/bin/env Rscript
source("R/data-process-angio.R")
saveRDS(build.data("seedMass"), "output/data-seedMass.rds")
