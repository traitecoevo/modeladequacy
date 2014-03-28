#!/usr/bin/env Rscript
source("R/data-process-angio.R")
saveRDS(build.data("sla"), "output/data-sla.rds")
