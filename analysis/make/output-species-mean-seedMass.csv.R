#!/usr/bin/env Rscript
## TODO: What is the slow step here?
library(methods)
suppressMessages(library(dplyr))
source("R/load-scrubbing-tools.R")

kew <- read.csv("data/kew.csv", stringsAsFactors=FALSE)

dat <- data.frame(gs=scrub.wrapper(kew$species),
                  seedMass=kew$value,
                  stringsAsFactors=FALSE)

errors <- get.errors()
kew.errors <- errors[errors$Dataset=="kewSeed" &
                     errors$trait=="seedMass",]
kew.errors$Original <- as.numeric(kew.errors$Original)
kew.errors$Changed  <- as.numeric(kew.errors$Changed)
kew.errors$genus_species <- gsub(" ", "_", kew.errors$genus_species)
#double column matching
replace.matrix <- which(
  outer(kew.errors$genus_species, dat$gs, "==") &
  outer(round(kew.errors$Original), round(dat$seedMass), "=="),
  arr.ind=TRUE)
## TODO: Again, only 21 of 60 match.
dat$seedMass[replace.matrix[,2]] <-
  kew.errors$Changed[replace.matrix[,1]]

#using modified plant list synonmy
pl.mod <- get.synonyms()
temp <- pl.mod$species[match(dat$gs, pl.mod$synonym)]
dat$gs[dat$gs %in% pl.mod$synonym] <- temp[!is.na(temp)]

## Build a little data frame with the species names and geometric mean
## of the trait:
dat.spp <-
  dat[complete.cases(dat),] %.%
  group_by(gs)              %.%
  summarise(mean = geometric.mean(seedMass),
            sd   = geometric.sd(seedMass))

write.csv(dat.spp, "output/species-mean-seedMass.csv", row.names=FALSE)

## TODO: Same issues as leafN.
sd <- mean(log10(dat.spp$sd[log10(dat.spp$sd) > 0.0001]), na.rm=TRUE)

## NOTE: Temporary checking section:
##
## Here I'm reading from data created by the old script, rather than
## the previous set of data; we agree with the generated set using the
## previous script, even if not with the previous data set.  I suspect
## new data are the cause there.
dat <- read.csv("output/species-mean-seedMass.csv",
                stringsAsFactors=FALSE)
ref <- read.csv("output/ref/species-mean-seedMass-oldscript.csv",
                stringsAsFactors=FALSE)
names(ref) <- c("gs", "mean")

library(testthat)
expect_that(dim(dat[1:2]), equals(dim(ref)))
expect_that(dat$gs,        equals(ref$gs))
expect_that(dat$mean,      equals(ref$mean))
expect_that(rownames(dat), equals(rownames(ref)))
