## internal functions for defining generic fxns and parameter names

modeltype <- function(fit, ...)
  UseMethod("modeltype")
modeldata <- function(fit, ...)
  UseMethod("modeldata")
modelpars <- function(fit, ...)
  UseMethod("modelpars")
modelinfo <- function(fit, ...)
  UseMethod("modelinfo")



modelphylo <- function(x, ...)
    UseMethod("modelphylo")



parnames.bm <- function()
  c("sigsq", "z0", "SE")
parnames.ou <- function()
  c("alpha", "sigsq", "z0", "SE")
parnames.eb <- function()
  c("a", "sigsq", "z0", "SE")
parnames.lambda <- function()
  c("lambda", "sigsq", "z0", "SE")
parnames.kappa <- function()
  c("kappa", "sigsq", "z0", "SE")
parnames.delta <- function()
  c("delta", "sigsq", "z0", "SE")
parnmaes.trend <- function()
  c("slope", "sigsq", "z0", "SE")
parnames.white <- function()
  c("sigsq", "z0", "SE")
