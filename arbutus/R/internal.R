## internal functions for defining generic fxns and parameter names

model.type <- function(fit, ...)
  UseMethod("model.type")
model.data <- function(fit, ...)
  UseMethod("model.data")
model.pars <- function(fit, ...)
  UseMethod("model.pars")
model.info <- function(fit, ...)
  UseMethod("model.info")






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
