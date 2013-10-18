## functions for pulling apart phylolm objects

model.type.phylolm <- function(fit, ...){
    model <- fit$model

    ## treat random and fixed root OU model as the same for now
    if (model == "OUrandomRoot" | model == "OUfixedRoot")
        model <- "OU"

    model
}


## It does not appear that the original data and phy is included
## Must include them separately for now...

model.data.phylolm <- function(fit, phy, data, ...){
    list(phy=phy, data=data)
}


model.pars.phylolm <- function(fit, ...){
    sigsq <- fit$sigma2

    ## match model to parameters
    model <- model.type(fit)

    if (model == "BM"){
        pars <- list(sigsq=sigsq, z0=NA, SE=0)
    } else {
        pars <- list(sigsq=sigsq, par=fit$optpar, z0=NA, SE=0)

        mm <- c("OU", "EB", "lambda", "kappa", "delta")
        par <- c("alpha", "a", "lambda", "kappa", "delta")
        tmp <- par[which(mm == model)]
        names(pars)[names(pars) == "par"] <- tmp
    }
    pars
}

#' @method model.info phylolm
#' @S3method model.info phylolm
model.info.phylolm <- function(fit, phy, data, ...){
    m <- list(data=model.data(fit, phy, data),
              pars=model.pars(fit),
              type=model.type(fit))
    class(m) <- "fitC"
    m
}
