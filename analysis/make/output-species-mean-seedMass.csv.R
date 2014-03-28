#!/usr/bin/env Rscript
library(methods)
suppressMessages(library(dplyr))
source("R/data-process-taxonomic.R")
source("R/data-process-angio.R")

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
dat <- update.synonomy(dat)

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
