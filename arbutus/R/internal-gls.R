
## fxns for parsing a gls object from fitting a pgls model using gls (from nlme)
## note here that the relevant data is the residuals, not the original trait values

model.type.gls <- function(fit, ...){

    ## figure out what model was used using ape's cor scheme
    ## convert model to geiger's nomenclature
    cor <- class(fit$modelStruct$corStruct)[1]

    if (is.null(cor) || !is.character(cor) || length(cor) != 1)
        stop("Failed to extract model type from gls fit")

    gmodels <- c("BM", "OU", "EB", "lambda")
    cormodels <- c("corBrownian", "corMartins", "corBlomberg", "corPagel")

    if (cor %in% cormodels == FALSE)
        stop("Could not identify correlation structure in gls fit")

    model <- gmodels[which(cormodels == cor)]

    model
}

model.data.gls <- function(fit, ...){
    phy <- attr(fit$modelStruct$corStruct, "tree")

    ## data is the residuals
    data <- as.numeric(resid(fit))
    names(data) <- phy$tip.label

    list(phy=phy, data=data)
}

model.pars.gls <- function(fit, ...){
    ## get model
    model <- model.type(fit)

    ## get phy
    phy <- model.data(fit)$phy

    ## pull out residuals
    rr <- resid(fit)
   

    if (model == "BM"){
       sigsq <- (rr %*% solve(vcv(phy)) %*% rr) / (Ntip(phy)-1) 
       pars <- list(sigsq=sigsq, z0=NA, SE=0)
    }
    
    ## if model is not BM, need to rescale the tree first, before computing sigsq
    if (model == "OU"){
        alpha <- as.numeric(fit$modelStruct$corStruct)
        tmp <- list(sigsq=1, alpha=alpha, SE=0)
        rescalephy <- model.phylo.ou(phy, pars=tmp)
        sigsq <- (rr %*% solve(vcv(rescalephy)) %*% rr) / (Ntip(phy)-1) 
        pars <- list(alpha=alpha, sigsq=sigsq, z0=NA, SE=0)
    }

    if (model == "EB"){
        a <- as.numeric(fit$modelStruct$corStruct)
        tmp <- list(sigsq=1, a=a, SE=0)
        rescalephy <- model.phylo.eb(phy, pars=tmp)
        sigsq <- (rr %*% solve(vcv(rescalephy)) %*% rr) / (Ntip(phy)-1) 
        pars <- list(a=a, sigsq=sigsq, z0=NA, SE=0)
    }

    if (model == "lambda"){
        lambda <- as.numeric(fit$modelStruct$corStruct)
        tmp <- list(sigsq=1, lambda=lambda, SE=0)
        rescalephy <- model.phylo.lambda(phy, pars=tmp)
        sigsq <- (rr %*% solve(vcv(rescalephy)) %*% rr) / (Ntip(phy)-1) 
        pars <- list(lambda=lambda, sigsq=sigsq, z0=NA, SE=0)
    }
    
    pars
}

model.info.gls <- function(fit, ...){
    m <- list(data=model.data(fit),
         pars=model.pars(fit),
         type=model.type(fit))
    class(m) <- "fitC"
    m
}
