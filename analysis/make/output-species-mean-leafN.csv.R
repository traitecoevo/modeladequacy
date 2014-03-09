#!/usr/bin/env Rscript
library(methods, quietly=TRUE)
suppressMessages(library(dplyr))
source("R/load-scrubbing-tools.R")

dat <- read.csv("data/wright-2004.csv", stringsAsFactors=FALSE)
dat <- dat[c("Species", "N.mass")]
names(dat) <- c("gs", "N.mass")

## Sort out synonomy and mispellings in species names.
dat$gs <- scrub.wrapper(dat$gs)

## Correct errors:
##
## [no known errors in glopnet, so nothing here]

## TODO: What is going on here?
## TODO: pl.mod should not be global
## TODO: wrap this into a function.
## using modified plant list synonmy
pl.mod <- get.synonyms()
temp <- pl.mod$species[match(dat$gs, pl.mod$synonym)]
dat$gs[dat$gs%in%pl.mod$synonym] <- temp[!is.na(temp)]

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

## NOTE: Temporary checking section:
dat <- read.csv("output/species-mean-leafN.csv",
                stringsAsFactors=FALSE)
ref <- read.csv("output/ref/species-mean-leafN.csv",
                stringsAsFactors=FALSE)
names(ref) <- c("gs", "mean")

library(testthat)
expect_that(dim(dat[1:2]), equals(dim(ref)))
expect_that(dat$gs,        equals(ref$gs))
expect_that(dat$mean,      equals(ref$mean))
expect_that(rownames(dat), equals(rownames(ref)))
