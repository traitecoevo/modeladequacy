#!/usr/bin/env Rscript
source("R/read-data-functions.R")
saveRDS(build.data("sla"), "output/data-sla.rds")
