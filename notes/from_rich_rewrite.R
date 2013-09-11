## The utility functions for parsing a fitContinuous fit would not get
## exported.
##
## As an aside, it looks like 'argn' is based off
## diversitree::argnames?  I do wish you guys had chatted with me
## before embarking on the rewrite; we could have factored a bunch of
## stuff into a support package and started something way more useful
## than the current set of duplication.
##
## NOTE: Given that geiger is not solidified yet, let's tweak
## fitContinuous to add a 'model' field to avoid the ugliness below.
parsefitContinuous <- function(fit){
  if (!inherits(fit, "gfit"))
    stop("Object does not look like a geiger fit")

  ## NOTE: How do we know that we only want the first name?  What
  ## about approaches that use more than one parameter?
  arg.names <- argn(fit$lik)[1]

  ## This section is entirely working out what model we have.  Would
  ## be easier if we could do
  ##   model <- fit$model
  all.models <- c("BM", "OU", "EB", "lambda", "kappa", "delta", "trend")
  param <- c("sigsq", "alpha", "a", "lambda", "kappa", "delta", "slope")
  tmp <- cbind(all.models, param)
  model <- as.character(tmp[tmp[,"param"] == arg.names, "all.models"])

  ## NOTE: This assumes that we're working with the BM model only,
  ## right?
  sigsq <- fit$opt$sigsq
  ## So, I've documented that requirement with an error condition for
  ## the time being:
  if (model != "BM")
    stop("Only BM models currently supported")
  
  z0 <- fit$opt$z0

  ## NOTE: Can you help me out?  What's going on here?  I would have
  ## throught that sigsq would be the model paramters (especially
  ## given that the name for the model parameters is going there).
  if (model == "BM"){ # BM does not have a model specific parameter
    model.p <- 0
  } else {
    model.p <- fit$opt[[1]]
  }

  ## NOTE: When is SE *not* included with the output?  And if it's not
  ## included is NA not a better value to take?
  ##
  ## NOTE: On second read through, this looks like it might be the
  ## error estimates for tips?  Is that right?
  if ("SE" %in% names(fit$opt)){
    SE <- fit$opt$SE
  } else {
    SE <- 0
  }

  ## NOTE: I think that ultimately a list with elements:
  ##   model -- character string (1)
  ##   par   -- parameters (>= 1)
  ##   SE    -- error estimates (length(par))
  ## will be more flexible
  par <- c(sigsq, z0, model.p, SE)
  names(par) <- c("sigsq", "z0", arg.names, "SE")

  list(model=model, par=par)
}

## Better still, geiger could provide something like this for pulling
## the coefficients out nicely:
coef.gfit <- function(object, with.se=FALSE, ...) {
  arg.names <- argn(object$lik)
  pars <- unlist(object$opt[arg.names])
  if (with.se) {
    attr(pars, "SE") <- object$opt$SE
  }
  pars
}

## This is really something that geiger should provide as an internal
## function -- then we can just assume that the interface exists.
fetchDataFitContinuous <- function(fit) {
  e <- environment(fit$lik)
  list(phy=get("phy", e), data=get("dat", e))
}

## So that we end up with:
parseFitContinuous <- function(fit) {
  if (!inherits(fit, "gfit"))
    stop("Object does not look like a geiger fit")
  # set as default for now.
  if (is.null(fit$model))
    fit$model <- "BM"

  model <- fit$model
  pars <- coef(fit, TRUE)
  se <- attr(pars, "SE")

  data <- fetchDataFitContinuous(fit)

  list(model=model,
       pars=c(pars), # the c() drops attributes
       se=se,
       # Augment with the underlying data, too
       phy=data$phy,
       data=data$data)
}

## Let's talk about this at some point -- this is where we end up
## thinking about per-branch functions, right?
rescale.fitContinuous <- function(obj) {
  if (obj$model == "BM")
    rescale.bm(obj$phy, obj$par[1], obj$SE)
  else
    stop("Not yet implemented")
}

## I've pulled this into its own function for now; we're likely to see
## this again?  Even if not, it will help keep the other function
## cleanish.
rescale.bm <- function(phy, sigsq, se) {
  unit.tree <- phy
  
  unit.tree$edge.length <- unit.tree$edge.length * sigsq

  if (!is.null(se)) {
    tips <- edge.is.terminal(unit.tree)
    unit.tree$edge.length[tips] <- unit.tree$edge.length[tips] + se^2
  }

  unit.tree
}
  
edge.is.terminal <- function(phy) {
  phy$edge[,2] <= Ntip(phy)
}
