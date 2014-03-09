#!/usr/bin/env Rscript
library(methods)
source("R/read-data-functions.R")

tree <- get.tree()

corrections <- read.delim("data/names-tr.txt", stringsAsFactors=FALSE)
corrections <- corrections[!(corrections$originalName %in%
                             sub("_"," ",tree$tip.label)),]

saveRDS(corrections, "output/corrections.rds")
