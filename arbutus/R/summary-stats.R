#' Calculate summary statistics on a unit tree
#'
#' @param unit.tree a 'unit.tree' object or list of unit.trees
#' @param stats a named list of summary statistics to calculate on the unit.tree.
#'   If no \code{stats} argument supplied, default summary statistics are used.
#'
#' @return A data.frame with the summary statistics.
#'
#' @note This function can be applied to either a single 'unit.tree' of object or a list of 'unit.tree' objects.
#'   If \code{stats=NULL} default summary statistics are used (see \code{\link[arbutus]{def.summ.stats} for details).
#'   User defined summary statistics can be supplied as a named list of functions (see examples). The functions supplied
#'   must take a unit.tree as argument and perform an operation on at least one of the elements of the object
#'   (see \code{link}{as.unit.tree})
#'  
#'
#' @seealso \code{\link[arbutus]{def.summ.stats}}, \code{\link{sigsq.reml}}, \code{\link{var.pic}}
#'   \code{\link{slope.con.bl}}, \code{\link{slope.con.asr}}, \code{\link{slope.con.nh}}
#'
#' @export summ.stats
#'
#' @examples
#'   data(geospiza)
#'   td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#'   phy <- td$phy
#'   dat <- td$data[,"wingL"]
#'   unit.tree <- as.unit.tree(phy, dat)
#'
#'   ## use default statistics
#'   summ.stat <- summ.stats(unit.tree, stats=NULL)
#'   summ.stat
#'
#'   ## user defined statistics
#'   mean.con <- function(x) mean(x$pics[,"contrasts"])
#'   max.con <- function(x) max(x$pics[,"contrasts"])
#'
#'   summ.stat.user <- summ.stats(unit.tree, stats=list(mean = mean.con, max = max.con))
#'   summ.stat.user
#'  
summ.stats <- function(unit.tree, stats=NULL){
  stats <- check.summ.stats(stats)
  if (inherits(unit.tree, "unit.tree")){
    # One unit.tree
    res <- do.call(cbind.data.frame,
                   lapply(stats, function(x) x(unit.tree)))
  } else if (is.list(unit.tree) && length(unit.tree) > 0) {
    # Multiple unit.trees
    res <- do.call(rbind, lapply(unit.tree, summ.stats, stats))
  } else {
    # Anything else
    stop("'unit.tree' must be a single unit.tree, or a list of them")
  }
  res
}

## arbutus:::sigsq.reml

## A function for computing the REML estimate of sigma^2
## Equal to the mean of the squared contrasts
## Takes a unit.tree object

sigsq.reml <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)

    mean(unit.tree$pics[,"contrasts"]^2)
}





## arbutus:::ks.pic

## A function for comparing the distribution of the contrasts
## to a normal distribution with mean 0 and sd = sqroot(mean(pics^2))
## uses D-statistic from KS test

ks.pic <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
        
    ## simulate a null distribution
    sd <- sqrt(mean(unit.tree$pics[,"contrasts"]^2))
    nd <- rnorm(10000, mean=0, sd=sd)

    ## KS test
    ksbm <- ks.test(unit.tree$pics[,"contrasts"], nd)$statistic

    unname(ksbm)
}







## arbutus:::var.pic

## A fxn for computing the variance in pic estimates
## Takes a unit.tree object

var.pic <- function(unit.tree) {
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
    
    var(abs(unit.tree$pics[,"contrasts"]))
}








## arbutus:::slope.pic.bl

## Fxn for computing the slope between the contrasts and their standard deviations
## Used as a measure of variation with respect to branch length

## Takes a unit.tree object

slope.pic.bl <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
    
    ## sd of the pics
    sd.pic <- sqrt(unit.tree$pics[,"variance"])

    ## absolute value of the pics
    abs.pic <- abs(unit.tree$pics[,"contrasts"])

    ## fit a linear model
    cv <- lm(abs.pic ~ sd.pic)

    unname(coef(cv)["sd.pic"])
}









## arbutus:::slope.pic.nh

## Fxn for computing the slope between the contrasts and the height of the node at which the contrast is computed
## Used as a measure of variation with respect to time

## Takes a unit.tree object

# NOTES:
##
## 1. Need to swap out branching.times() with a function that works for non-ultrametric trees (working on the fxn nodeAge to do this)

slope.pic.nh <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)

    ## node heights (see note 1)
    nh <- branching.times(unit.tree$phy)

    ## absolute value of the pics
    abs.pic <- abs(unit.tree$pics[,"contrasts"])

    ## fit a linear model

    ch <- lm(abs.pic ~ nh)

    unname(coef(ch)["nh"])
}








## arbutus:::slope.pic.asr

## Fxn for computing the slope between the contrasts and the inferred ancestral state at the node
## Used as a measure of variation with respect to trait value

## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data



## NOTES:
##
## 1. Use "pic" method in ape. Not sure if this is appropriate or not.


slope.pic.asr <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
    
    ## reconstruct ancestors using pic method
    asr <- ace(unit.tree$data, unit.tree$phy, method="pic")$ace

    ## absolute value of pics
    abs.pic <- abs(unit.tree$pics[,"contrasts"])

    ## fit linear model
    ca <- lm(abs.pic ~ asr)

    unname(coef(ca)["asr"])
}









## arbutus:::def.summ.stats

## function to build list of default summary statistics
## used internally in traitStats

## takes no arguments

def.summ.stats <- function()
    list("reml.sigsq"=sigsq.reml, "var.con"=var.pic, "slope.con.var"=slope.pic.bl, "slope.con.asr"=slope.pic.asr, "slope.con.nh"=slope.pic.nh, "ks.dstat"=ks.pic)

check.summ.stats <- function(stats) {
  if (is.null(stats))
    stats <- def.summ.stats()
  else if (!is.list(stats) || is.null(names(stats)) ||
           any(names(stats) == "" || length(stats) == 0))
    ## (note that the second check does not always enforce
    ## length(stats) > 0)
    stop("'stats' must be a named list with at least one element")
  else if (!all(sapply(stats, is.function)))
    stop("All elements of stats must be functions")
  stats
}
