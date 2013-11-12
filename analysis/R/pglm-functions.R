library(ape)
library(diversitree)


default.argnames.pglm.bm <- function()
    c("y0", "slope")

default.argnames.pglm.lambda <- function()
    c("y0", "slope", "lambda")


## make info functions
make.info.pglm.bm <- function(phy){
    list(name="pglm.bm", name.pretty="PGLM", np = 2L,
         argnames=default.argnames.pglm.bm(), ny=3L, k=NA, idx.e=NA,
         idx.d=NA,  phy=phy, ml.default = "subplex", mcmc.lowerzero=TRUE,
         doc=NULL,
         reference=c("Freckleton R.P. 2012 Methods in Ecology and Evolution 3:940-947"))
}

make.info.pglm.lambda <- function(phy){
    list(name="pglm.lambda", name.pretty="PGLM with lambda", np = 3L,
         argnames=default.argnames.pglm.lambda(), ny=3L, k=NA, idx.e=NA,
         idx.d=NA,  phy=phy, ml.default = "subplex", mcmc.lowerzero=TRUE,
         doc=NULL,
         reference=c("Freckleton R.P. 2012 Methods in Ecology and Evolution 3:940-947"))
}


## check control
check.control.pglm <- function(control){
    defaults <- list(method="ML")
    if (length(control$method) != 1) 
        stop("control$method must be a scalar")

    methods <- c("ML", "REML")
    if (!(control$method %in% methods)) 
        stop(sprintf("control$method must be in %s", paste(methods, 
            collapse = ", ")))

    control
}





## make cache
make.cache.pglm.bm <- function(tree, states.x, states.y, control){
    method <- control$method
    tree <- diversitree:::check.tree(tree, ultrametric=FALSE)
    cache <- diversitree:::make.cache(tree)
    cache$states.x <- diversitree:::check.states(tree, states.x, as.integer=FALSE)
    cache$states.y <- diversitree:::check.states(tree, states.y, as.integer=FALSE)
    cache$pic.x <- pic(x=cache$states.x, phy=tree, var.contrasts = TRUE)
    cache$pic.y <- pic(x=cache$states.y, phy=tree, var.contrasts = TRUE)
    cache$info <- make.info.pglm.bm(tree)
    cache
}
    


## make log likelihood eqn for ML estimate






make.pglm.bm <- function(tree, states.x, states.y, control=list()){
    control <- check.control.pglm(control)
    cache <- make.cache.pglm.bm(tree, states.x, states.y, control)

    if (control$method == "ML")
        all.branches <- make.all.branches.pglm.ml(cache, control)

    if (control$method == "REML")
        all.branches <- make.all.branches.pglm.reml(cache, control)

    ll <- function(pars){
        check.pars.pglm.bm(pars)
        ans <- all.branches(pars)
    }
    class(ll) <- c("pglm.bm", "dtlik", "function")
    ll
}
    

