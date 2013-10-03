## Arrange to quietly load these packages so that messages don't
## overwhelm the output; lots of dependent packages and some conflicts
## generate a lot of noise.
library(testthat)
suppressMessages(library(arbutus))
suppressMessages(library(diversitree))
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
