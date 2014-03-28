#!/usr/bin/env Rscript
library(methods, quietly=TRUE)
suppressMessages(library(dplyr))
source("R/data-process-taxonomic.R")
source("R/data-process-angio.R")

dat <- read.csv("data/wright-2004.csv", stringsAsFactors=FALSE)
dat <- dat[c("Species", "N.mass")]
names(dat) <- c("gs", "N.mass")

## Sort out synonomy and mispellings in species names.
dat$gs <- scrub.wrapper(dat$gs)

## Correct errors:
##
## [no known errors in glopnet, so nothing here]

dat <- update.synonomy(dat)

## Build a little data frame with the species names and geometric mean
## of the trait:
dat.spp <-
  dat[complete.cases(dat),] %.%
  group_by(gs)              %.%
  summarise(mean = geometric.mean(N.mass),
            sd   = geometric.sd(N.mass))

write.csv(dat.spp, "output/species-mean-leafN.csv", row.names=FALSE)

## TODO: Save this value to a file somewhere, and document more fully
## what it is.  It's a bit confusing because at this point we have raw
## values for leafN but we are going to eventually work in log-10
## space.  I've possibly sorted that out longer term by writing out
## the sd values with the data above.
sd <- mean(log10(dat.spp$sd[log10(dat.spp$sd) > 0.0001]), na.rm=TRUE)
