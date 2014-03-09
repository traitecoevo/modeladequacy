#!/usr/bin/env Rscript
library(methods)
suppressMessages(library(dplyr))
source("R/load-scrubbing-tools.R")

## Two data sets
wright <- read.csv("data/wright-2004.csv",
                   stringsAsFactors=FALSE)
# I am at a loss why this is coming through as character.
wright$LogLMA <- as.numeric(wright$LogLMA)

# Extract just the columns we want, converting lma to sla
wright <- data.frame(gs=wright$Species,
                     sla=10000/10^(wright$LogLMA),
                     dataset="glop",
                     stringsAsFactors=FALSE)

leda <- read.csv("data/leda.csv",
                 stringsAsFactors=FALSE)

leda <- data.frame(gs=leda$SBS.name,
                   sla=10 * leda$SLA.mean,
                   dataset="leda",
                   stringsAsFactors=FALSE)

dat <- rbind(wright, leda)

## Sort out synonomy and mispellings in species names.
dat$gs <- scrub.wrapper(dat$gs)

## Correct errors:
##
errors <- get.errors()
leda.errors <- errors[errors$Dataset=="LEDA" & errors$trait=="sla",]
leda.errors$Original <- as.numeric(leda.errors$Original)
leda.errors$Changed  <- as.numeric(leda.errors$Changed)
leda.errors$gs <- sub(" ", "_", leda.errors$genus_species, fixed=TRUE)

# double column matching

## Hmm, looks like this was not working.  I needed to do the
## underscore substitute, but even then I don't get all the matches I
## need.  I'd expect 21 here, but only get 10.
##
## TODO: talk through this with Matt, replace with something more
## robust?
replace.matrix <- which(
  outer(leda.errors$gs, dat$gs, "==") &
  outer(round(leda.errors$Original), round(dat$sla), "=="),
  arr.ind=TRUE)
dat$sla[replace.matrix[,2]] <-
  leda.errors$Changed[replace.matrix[,1]]

## TODO: See mean-leafN
## using modified plant list synonmy
pl.mod <- get.synonyms()
temp <- pl.mod$species[match(dat$gs, pl.mod$synonym)]
dat$gs[dat$gs %in% pl.mod$synonym] <- temp[!is.na(temp)]

## Build a little data frame with the species names and geometric mean
## of the trait:
dat.spp <-
  dat[complete.cases(dat),] %.%
  group_by(gs)              %.%
  summarise(mean = geometric.mean(sla),
            sd   = geometric.sd(sla))

write.csv(dat.spp, "output/species-mean-sla.csv", row.names=FALSE)

## TODO: Same issues as leafN.  However, because I could not run the
## SLA code, this value is unchecked.
sd <- mean(log10(dat.spp$sd[log10(dat.spp$sd) > 0.0001]), na.rm=TRUE)

## NOTE: Temporary checking section:
##
## TODO: This differs very very slightly.  We've lost two species, but
## the numbers largely look right.  Look into this later once I've
## spoken with Matt, but this really requires that the original
## version can run.
dat <- read.csv("output/species-mean-sla.csv",
                stringsAsFactors=FALSE)
ref <- read.csv("output/ref/species-mean-sla.csv",
                stringsAsFactors=FALSE)
names(ref) <- c("gs", "mean")

if (FALSE) {
  library(testthat)
  expect_that(dim(dat[1:2]), equals(dim(ref)))
  expect_that(dat$gs,        equals(ref$gs))
  expect_that(dat$mean,      equals(ref$mean))
  expect_that(rownames(dat), equals(rownames(ref)))
}
