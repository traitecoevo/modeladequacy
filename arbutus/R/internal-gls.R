
## fxns for parsing a gls object from fitting a pgls model using gls (from nlme)
## note here that the relevant data is the residuals, not the original trait values

modeltype.gls <- function(fit, ...){

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

modeldata.gls <- function(fit, ...){
    phy <- attr(fit$modelStruct$corStruct, "tree")

    ## data is the residuals
    data <- as.numeric(resid(fit))
    names(data) <- phy$tip.label

    list(phy=phy, data=data)
}

modelpars.gls <- function(fit, ...){
    ## get model
    model <- modeltype(fit)

    ## get sigsq
    phy <- modeldata(fit)$phy
    rr <- resid(fit)
    sigsq = (rr %*% solve(vcv(phy)) %*% rr) / (Ntip(phy)-1)

    if (model == "BM")
       pars <- list(sigsq=sigsq, SE=0)

    if (model == "OU")
        pars <- list(alpha=as.numeric(fit$modelStruct$corStruct), sigsq=sigsq, , SE=0)

    if (model == "EB")
        pars <- list(a=as.numeric(fit$modelStruct$corStruct), sigsq=sigsq, SE=0)

    if (model == "lambda")
        pars <- list(lambda=as.numeric(fit$modelStruct$corStruct),sigsq=sigsq, SE=0)

    pars
}

modelinfo.gls <- function(fit, ...){
    m <- list(data=modeldata(fit),
         pars=modelpars(fit),
         type=modeltype(fit))
    class(m) <- "model.fitC"
    m
}
