## fxns for parsing a caper 'pgls' object
## note here that the relevant data is the residuals,
## not the original trait values

## this is awkward at the moment
## if parameter not estimated, still returned as equal to 1
## but if it is estimated to be 1, save output returned
## so hacked this together
model.type.pgls <- function(fit, ...){
    types <- c("kappa", "lambda", "delta")
    tmp <- fit$param.CI
    chk <- Filter(Negate(is.null), tmp)
    are.fit <- types %in% names(chk)

    ## if no parameters are estimated, BM
    if (all(are.fit == FALSE)){
        model <- "BM"
    } else {
        model <- types[are.fit]
        
        ## if multiple models estimated, return error
        if (length(model) > 1)
            stop("multiple models fit. cannot assess adequacy")

    }
    model
}


model.data.pgls <- function(fit, ...){
    phy <- fit$data$phy

    ## data is the residuals
    data <- as.numeric(resid(fit))
    names(data) <- phy$tip.label

    list(phy=phy, data=data)
}


## See issues #50 and #51
## Also note that REML estimators are not available for this pkg
model.pars.pgls <- function(fit, ...){
    model <- model.type(fit)
    if (model == "BM"){
        p <- NULL
    } else {
        p <- fit$param[model]
    }
    pars <- as.list(c(p, sigsq=estimate.sigma2.pgls(fit), SE=0, z0=NA))
}



estimate.sigma2.pgls <- function(fit, ...){
    phy <- model.data(fit)$phy

    pars <- get.pgls.pars.unity(fit)
    phy <- model.phylo.rescale(model.type(fit))(phy, pars)

    ## The data that we care about are the residuals of the model fit;
    ## it is these that are assumed to be distributed according to the
    ## tree's VCV.
  
    rr <- as.numeric(resid(fit))
    names(rr) <- phy$tip.label
    cmp <- resid(fit)
  
    phy.u <- as.unit.tree(phy, data=rr)
    ## This is what sigsq.est() (in summary-stats.R) is doing, but I
    ## (RGF) find it easier to think about if it's explicitly copied
    ## here.
    s2 <- mean(phy.u$pics[,"contrasts"]^2) # REML estimate

    ## now convert to ml
    ## NOTE: In contrast with gls, caper::pgls always fits ML rather
    ## than REML.
    s2 <- reml.to.ml(s2, Ntip(phy), length(coef(fit)))

    s2
}


## helper function for pulling out parameters for estimate sigma2.pgls
get.pgls.pars.unity <- function(fit){
    model <- model.type(fit)
    if (model == "BM"){
        pars <- c(list(sigsq=1), list(SE=0))
    } else {
        pars <- c(list(sigsq=1), as.list(fit$param[model]), list(SE=0))
    }
    pars
}


#' @method model.info pgls
#' @S3method model.info pgls
model.info.pgls <- function(fit, ...){
    m <- list(data=model.data(fit),
              pars=model.pars(fit),
              type=model.type(fit))
    class(m) <- "fitC"
    m
}
    
    
        
    
