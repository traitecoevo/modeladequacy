#!/usr/bin/env Rscript
source("R/data-process-angio.R")
saveRDS(build.data("leafN"), "output/data-leafN.rds")
