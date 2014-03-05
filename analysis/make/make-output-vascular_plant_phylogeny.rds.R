#!/usr/bin/env Rscript
library(ape)
saveRDS(read.tree("data/vascular_plant_phylogeny.tre"),
        "output/vascular_plant_phylogeny.rds")
