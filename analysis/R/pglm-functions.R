library(ape)
library(diversitree)


default.argnames.pglm.bm <- function()
    c("slope", "epsilon")

default.argnames.pglm.lambda <- function()
    c("slope", "epsilon", "lambda")


## make info functions
make.info.pglm.bm <- function(phy){
    list(name="pglm.bm", name.pretty="PGLM", np = 1L,
         argnames=default.argnames.pglm.bm(), ny=3L, k=NA, idx.e=NA,
         idx.d=NA,  phy=phy, ml.default = "subplex", mcmc.lowerzero=TRUE,
         doc=NULL,
         reference=c("Freckleton R.P. 2012 Methods in Ecology and Evolution 3:940-947"))
}

make.info.pglm.lambda <- function(phy){
    list(name="pglm.lambda", name.pretty="PGLM with lambda", np = 2L,
         argnames=default.argnames.pglm.lambda(), ny=3L, k=NA, idx.e=NA,
         idx.d=NA,  phy=phy, ml.default = "subplex", mcmc.lowerzero=TRUE,
         doc=NULL,
         reference=c("Freckleton R.P. 2012 Methods in Ecology and Evolution 3:940-947"))
}


## check control
check.control.pglm <- function(control){
    defaults <- list(method="ML")
    control <- modifyList(defaults, control)
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



## Frecklton's pic function
## NOT SURE EXACTLY WHY HE USES THIS
pic.vcv <- function (x, phy, scaled = TRUE, var.contrasts = TRUE) 
{
    if (class(phy) != "phylo") 
        stop("object 'phy' is not of class \"phylo\"")
        
    if (is.null(phy$edge.length)) 
        stop("your tree has no branch lengths: you may consider setting them equal to one, or using the function `compute.brlen'.")
       
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    
    if (nb.node != nb.tip - 1) 
        stop("'phy' is not rooted and fully dichotomous")
        
    if (length(x) != nb.tip) 
        stop("length of phenotypic and of phylogenetic data do not match")
        
    if (any(is.na(x))) 
        stop("the present method cannot (yet) be used directly with missing data: you may consider removing the species with missing data from your tree with the function `drop.tip'.")
        
    phy <- reorder(phy, "pruningwise")
    
    phenotype <- numeric(nb.tip + nb.node)
    
    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    }
    else {
        if (all(names(x) %in% phy$tip.label)) 
            phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning("the names of argument \"x\" and the names of the tip labels did not match: the former were ignored in the analysis.")
        }
    }
    
    contr <- var.con <- numeric(nb.node)
    
    ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node), 
        as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]), 
        as.double(phy$edge.length), as.double(phenotype), as.double(contr), 
        as.double(var.con), as.integer(var.contrasts), as.integer(scaled), 
        PACKAGE = "ape")
        
    contr <- ans[[7]]
    
    if (var.contrasts) {
        contr <- cbind(contr, ans[[8]])
        dimnames(contr) <- list(1:nb.node + nb.tip, c("contrasts", 
            "variance"))
    } else names(contr) <- 1:nb.node + nb.tip
        
   	idx <- which(ans[[3]] == (nb.tip+1) )
    root.v <- ans[[5]][idx]
    V <- root.v[1]*root.v[2]/(sum(root.v))
    
    root.state <- ans[[6]][nb.tip + 1]
    
    return(list(contr = contr, root.v = root.v, V = V, ans = ans, root.state = root.state))
}




## function for calculating phylogenetic residuals
calc.phylo.resid <- function(u.x, u.y, b, vc, e){
    (u.y - u.x * b)^2 / (vc * e)
}




                             
## make log likelihood eqn for ML estimate
make.all.branches.pglm.bm.ml <- function(cache){
    u.x <- cache$pic.x[,"contrasts"]
    u.y <- cache$pic.y[,"contrasts"]
    vc <- cache$pic.x[,"variance"]
    n <- length(cache$info$phy$tip.label)

    ## likelihood function
    ll <- function(pars){
        b <- pars[1]
        e <- pars[2]
        -0.5 * (n * log(2*pi) + sum(log(vc)) +
                sum(sapply(seq_len(n-1), function(x)
                           calc.phylo.resid(u.x[x], u.y[x], vc[x], b, e))))
    }
    ll
}
    


check.pars.pglm.bm <- function(pars){
    if (length(pars) != 2)
        stop("Incorrect paramter length")
}




make.pglm.bm <- function(tree, states.x, states.y, control=list()){
    control <- check.control.pglm(control)
    cache <- make.cache.pglm.bm(tree, states.x, states.y, control)

    if (control$method == "ML")
        all.branches <- make.all.branches.pglm.bm.ml(cache)

    if (control$method == "REML")
        all.branches <- make.all.branches.pglm.bm.reml(cache)

    ll <- function(pars){
        check.pars.pglm.bm(pars)
        ans <- all.branches(pars)
        ans
    }
    class(ll) <- c("pglm.bm", "dtlik", "function")
    ll
}
    

