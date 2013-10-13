## internal functions for defining generic fxns and parameter names

model.type <- function(fit, ...)
  UseMethod("model.type")
model.data <- function(fit, ...)
  UseMethod("model.data")
model.pars <- function(fit, ...)
  UseMethod("model.pars")







#' @title Get model and data from model fitted object
#'
#' @description Generic function for pulling apart model fitted objects
#' and obtaining phylogeny, data, model type and parameters
#'
#' @param fit a model fitted object (see \code{\link{as.unit.tree}} for details)
#' @param ... additional arguments to be passed to \code{model.info}
#'
#' @details These functions are primarily used internally within the function
#' \code{\link{as.unit.tree}}. A different \code{model.info.xx} function is built for every implemented
#' class of model fitted object and is tailored to the specifics of the form of model output.
#' The return value will differ depending on the model used. For example, if a single phylogeny
#' is supplied which has been used to fit a single rate model, \code{model.info} will return
#' a object of class \code{fitC}, consisting of the original tree to be rescaled, the
#' original comparative dataset, the model class (e.g. "BM", "OU", etc.) and the model-specific
#' parameters.
#'
#' The returned object needs to be of a class for which a \code{\link{make.model.phylo}} function
#' exists. \code{\link{make.model.phylo}} will take the object output by \code{model.info} and
#' return a 'unit.tree' object.
#'
#' @return a classed object compatibile with a \code{\link{make.model.phylo}} function.
#'
#' @seealso \code{\link{as.unit.tree}}, \code{\link{make.model.phylo}}
#'
#' @export model.info
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' \dontrun{
#' require(geiger)
#' ## fit Brownian motion model
#' ## using geiger's fitContinuous function
#' fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
#'                                  control=list(niter=10))
#'
#' ## get model info using internal arbutus function
#' info.bm <- model.info(fit.bm)
#'
#' ## now using diversitree function
#' require(diversitree)
#' bmlik <- make.bm(phy, data)
#' fit.bm.dt <- find.mle(bmlik, x.init=1)
#'
#' info.bm.dt <- model.info(fit.bm.dt, lik=bmlik)
#' }
#' 
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
