#!/usr/bin/env Rscript
source("R/read-data-functions.R")
saveRDS(build.data("leafN"), "output/data-leafN.rds")
