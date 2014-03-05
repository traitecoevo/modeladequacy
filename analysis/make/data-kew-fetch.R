#!/usr/bin/env Rscript
source("R/kew.R")
dir.create("data/kew", FALSE)
kew.fetch.clades()
