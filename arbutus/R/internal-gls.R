
## fxns for parsing a gls object from fitting a pgls model using gls (from nlme)
## note here that the relevant data is the residuals, not the original trait values

model.type.gls <- function(fit, ...){
    cor <- class(fit$modelStruct$corStruct)[1]
    if (is.null(cor) || !is.character(cor) || length(cor) != 1)
        stop("Failed to extract model type from gls fit")

    ## Translate from ape's cor scheme convert model to geiger's:
    tr <- c(corBrownian = "BM",
            corMartins  = "OU",
            corBlomberg = "EB",
            corPagel    = "lambda")

    model <- tr[[cor]]
    if (is.na(model))
      stop("Could not identify correlation structure in gls fit")

    model
}

model.data.gls <- function(fit, ...){
    phy <- attr(fit$modelStruct$corStruct, "tree")

    ## data is the residuals
    data <- as.numeric(resid(fit))
    names(data) <- phy$tip.label

    list(phy=phy, data=data)
}

## See issues #50 and #51 here.
model.pars.gls <- function(fit, ...){
  p <- translate.ape.arbutus(coef(fit$modelStruct$corStruct, FALSE))
  as.list(c(p, sigsq=estimate.sigma2.gls(fit), SE=0, z0=NA))
}

#' @method model.info gls
#' @S3method model.info gls
model.info.gls <- function(fit, ...){
    m <- list(data=model.data(fit),
         pars=model.pars(fit),
         type=model.type(fit))
    class(m) <- "fitC"
    m
}

## TODO (RGF): I am not certain that EB is correct here, because the
## data that I fit with the EB model did not converge.  Should try
## this more carefully soon.
translate.ape.arbutus <- function(x) {
  tr <- c(a="g")
  idx <- match(names(x), tr)
  i <- !is.na(idx)
  names(x)[i] <- names(tr)[idx][i]
  x
}

## Helper function to compute sigma2 from a gls fit.  This will return
## either the REML or the ML sigma2, depending on what was originally
## chosen with the 'method' argument to gls().  The parameters are
## extracted from the correlation structure, though the parameters
## need a little tweaking to get.
##
## Note: I (MWP) changed rescale variable name to tree.rescale to
## avoid confusion with similar rescale.phylo fxn in geiger
estimate.sigma2.gls <- function(fit) {
  phy <- model.data(fit)$phy
  pars <- model.pars(fit$modelStruct$corStruct)

  ## Appropriately rescale the tree:
  tree.rescale <- switch(model.type(fit),
                    BM=model.phylo.bm,
                    OU=model.phylo.ou,
                    EB=model.phylo.eb,
                    lambda=model.phylo.lambda,
                    stop("Unknown model type ", model.type(fit)))
  phy <- tree.rescale(phy, pars)

  ## The data that we care about are the residuals of the model fit;
  ## it is these that are assumed to be distributed according to the
  ## tree's VCV.
  
  rr <- as.numeric(resid(fit))
  names(rr) <- names(resid(fit))
  cmp <- resid(fit)
  
  phy.u <- as.unit.tree(phy, data=rr)
  ## This is what sigsq.est() (in summary-stats.R) is doing, but I
  ## (RGF) find it easier to think about if it's explicitly copied
  ## here.
  s2 <- mean(phy.u$pics[,"contrasts"]^2) # REML estimate

  if (fit$method == "ML")
    s2 <- reml.to.ml(s2, fit$dims$N, fit$dims$p)

  s2
}


model.pars.corStruct <- function(fit, ...) {
  ## RGF: What happens when SE != 0; is that just not possible with
  ## this approach?
  pars <- c(list(sigsq=1), as.list(coef(fit, FALSE)), list(SE=0))
  translate.ape.arbutus(pars)
}

## TODO: Not completely sure about this, actually.  I think that we
## are actually looking for (n - (k - 1)).  See
## inst/tests/internal-gls-diversitree.R.
reml.to.ml <- function(s2, n, k)
  s2 * (n - k) / n  
ml.to.reml <- function(s2, n, k)
  s2 * n / (n - k)

