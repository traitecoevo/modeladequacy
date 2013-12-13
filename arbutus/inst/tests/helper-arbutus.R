## Arrange to quietly load these packages so that messages don't
## overwhelm the output; lots of dependent packages and some conflicts
## generate a lot of noise.
library(testthat)
suppressMessages(library(arbutus))
suppressMessages(library(geiger))
## note: geiger must come before diversitree due to conflicts
## between diversitree and coda
suppressMessages(library(diversitree))
library(nlme) # for gls
library(phylolm)
library(caper) 

library(parallel)

## Additional tests
is_greater_than <- function(value) {
  function(actual)
    expectation(actual > value, paste("is not greater than", value))
}

is_less_than <- function(value) {
  function(actual)
    expectation(actual < value, paste("is not less than", value))
}

fitContinuousQuiet <- function(...)
  suppressWarnings(fitContinuous(...))

## Work around - I don't see how to do this with geiger.
rescale.se <- function(phy, ..., SE=0) {
  phy <- rescale(phy, ...)
  tips <- phy$edge[,2] <= Ntip(phy)
  phy$edge.length[tips] <- phy$edge.length[tips] + SE
  phy
}
