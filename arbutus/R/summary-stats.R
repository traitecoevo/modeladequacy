#' @title Calculate summary statistics on a unit tree
#'
#' @description Calculates a set of summary statistics on the contrasts included with a 'unit.tree' object.
#'
#' @param unit.tree a 'unit.tree' object or list of unit.trees
#' @param stats a named list of summary statistics to calculate on the unit.tree.
#'   If no \code{stats} argument supplied, default summary statistics are used.
#'
#' @details This function can be applied to either a single 'unit.tree' of object or a list of 'unit.tree' objects.
#' If \code{stats=NULL} default summary statistics are used. The default summary statistics are the following:
#'  \enumerate{
#'   \item{reml.sigsq: }{The mean of the squared contrasts. This is equivalent to the REML estimate of sigsq.}
#'
#'   \item{var.con: }{The variance of the absolute value of the contrasts.}
#'
#'   \item{slope.con.var: }{The slope of a linear model fit between the contrasts and their expected variances.}
#'
#'   \item{slope.con.asr: }{The slope of a linear model fit between the contrasts and their inferred ancestral state.}
#'
#'   \item{slope.con.nh: }{The slope of a linear model fit between the contrasts and the node height at which they were calculated.}
#'
#'   \item{ks.pic: }{The D-statistic from a KS test comparing the distribution of the contrasts to a normal distribution  with mean 0 and variance equal to the square root of the squared mean of the contrasts.}
#'  }
#' User defined summary statistics can be supplied as a named list of functions (see examples). The functions supplied
#' must take a unit.tree as argument and perform an operation on at least one of the elements of the object
#' (see \code{\link{as.unit.tree}} for details).
#'
#' @return A data.frame with the calculated summary statistics across all unit.trees provided.
#'
#'
#' @export summ.stats
#'
#' @seealso \code{\link{def.summ.stats}}, \code{\link{sigsq.reml}}, \code{\link{var.pic}}, \code{\link{slope.con.bl}}, \code{\link{slope.con.nh}}, \code{\link{slope.con.asr}}, \code{\link{ks.pic}}
#'
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## use default statistics
#' summ.stat <- summ.stats(unit.tree, stats=NULL)
#' summ.stat
#'
#' ## user defined statistics
#' mean.con <- function(x) mean(x$pics[,"contrasts"])
#' max.con <- function(x) max(x$pics[,"contrasts"])
#'
#' summ.stat.user <- summ.stats(unit.tree,
#'                     stats=list(mean = mean.con, max = max.con))
#' summ.stat.user
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









#' @title REML estimate of BM parameter
#'
#' @description Calculates the mean of the squared contrasts,
#' which is equal to the REML estimate of the BM rate parameter.
#' Used as a summary statistic to evaluate model adequacy
#'
#' @param unit.tree a 'unit.tree' object
#'
#' @details This summary statistic is used to assess the ability of the model
#' to adequate estimate the rate of evolution. On datasets simulated on a unit.tree (see \code{\link{sim.char.unit}}),
#' the expected value is equal to 1.
#'
#' The REML estimate of \eqn{sigma^2} is included as a default
#' summary statistic in the function \code{\link{summ.stats}}
#'
#' @return estimate of \eqn{sigma^2}
#'
#' @export sigsq.reml
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{def.summ.stats}}
#'
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## Estimate sigsq from unit.tree
#'
#' sigsq.reml(unit.tree)
#' 
sigsq.reml <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)

    mean(unit.tree$pics[,"contrasts"]^2)
}








#' @title Evaluate normality of contrasts
#'
#' @description Performs Kolmogorov-Smirnov test comparing distribution
#' of contrasts with that of a normal distribution with mean 0 and variance
#' equal to the square root of the mean of squared contrasts
#'
#' @param unit.tree a unit.tree object
#'
#' @details This summary statistics is used to evaluate whether the assumption of multivariate normailty
#' is appropriate. If the model which generated the data is the fitted model, we expect the
#' square root ofthe mean of squared contrasts to be equal to 1. The empirical estimate is used
#' rather than assume a variance of 1 to reduce the overlap between this summary statistic and the
#' REML estimate of \eqn{sigma^2} (see \code{\link{reml.sigsq}}). The Kolmogorov-Smirnov (KS) test is a
#' non-parameteric test which computes the maximum distance \eqn{D} between two cumulative distribution functions.
#' Running the test multiple times on the same data will produce slightly different values due to the fact
#' that the null distribution is produced by randomly drawing from a normal distribution.
#'
#' The KS-D statistic is included as a default summary statistic
#' in the function \code{\link{summ.stats}}.
#'
#' @return the D-statistic from a KS-test
#'
#' @export ks.pic
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{def.summ.stats}}, \code{\link{stats::ks.test}}
#'
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## KS-D statistic
#' ks.pic(unit.tree)
#'
ks.pic <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
        
    ## simulate a null distribution
    sd <- sqrt(mean(unit.tree$pics[,"contrasts"]^2))
    nd <- rnorm(10000, mean=0, sd=sd)

    ## KS test
    ksbm <- suppressWarnings(ks.test(unit.tree$pics[,"contrasts"], nd)$statistic)

    unname(ksbm)
}







#' @title Variance of contrasts
#'
#' @description Calculates the variance of the absolute value of the contrasts
#'
#' @param unit.tree a 'unit.tree' object
#'
#' @details This summary is used to evaluate whether the model is sufficiently
#' capturing the variation in rates across the tree.
#'
#' The variance of contrasts is included as a default summary statistic
#' in the function \code{\link{summ.stats}}
#'
#' @return the estimate of the variance of the absolute value of the contrasts
#'
#' @export var.pic
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{def.summ.stats}}
#'
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## estimate variance of contrasts
#' var.pic(unit.tree)
#' 
var.pic <- function(unit.tree) {
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
    
    var(abs(unit.tree$pics[,"contrasts"]))
}








#' @title Relationship between contrasts and their variances
#'
#' @description Fits a linear model between the absolute value of constrasts
#' and their expected variances
#'
#' @param unit.tree a unit.tree object
#'
#' @details This summary statistic is the estimate of the slope from fitting
#' a linear model between the contrasts and their expected variance. It is used
#' to evaluate whether the model is adequately capturing variation relative to branch
#' lengths. If the generating model is correct, we expect there not to be a relationship
#' between the contrasts and their variances (i.e. slope ~ 0).
#'
#' The slope of the contrasts and their variances is included as a default summary statistic
#' in the function \code{\link{summ.stats}}
#'
#' @return the estimated slope paramter
#'
#' @export slope.pic.var
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{def.summ.stats}}, \code{\link{stats::lm}}
#'
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## estimate slope
#' slope.pic.var(unit.tree)
#' 
slope.pic.var <- function(unit.tree){
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









#' @title Relationship between contrasts and node height
#'
#' @description Fits a linear model between the absolute value of the contrasts
#' and the node height at which they were calculated
#'
#' @param unit.tree a 'unit.tree' object
#'
#' @details This summary statistic is the estimate of the slope from fitting
#' a linear model between the contrasts and the node height (i.e. tree depth)
#' at which they were calculated. It is used
#' to evaluate whether the model is adequately capturing variation relative to tree depth.
#' If the generating model is correct, we expect there not to be a relationship
#' between the contrasts and their node height (i.e. slope ~ 0).
#'
#' The slope of the contrasts and their variances is included as a default summary statistic
#' in the function \code{\link{summ.stats}} 
#' 
#' @return the estimated slope parameter
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{def.summ.stats}}, \code{\link{stats::lm}}
#'
#' @export slope.pic.nh
#'
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## estimate slope
#' slope.pic.nh(unit.tree)
#' 
slope.pic.nh <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)

    ## node heights 
    ht <- edge.height(unit.tree$phy)
    N <- Ntip(unit.tree$phy)
    nh <- ht$end[c((N + 1):nrow(ht))]
    names(nh) <- rownames(ht[c((N + 1):nrow(ht)),])

    ## absolute value of the pics
    abs.pic <- abs(unit.tree$pics[,"contrasts"])

    ## fit a linear model

    ch <- lm(abs.pic ~ nh)

    unname(coef(ch)["nh"])
}








#' @title Relationship between contrasts and inferred ancestral state
#'
#' @description Fits a linear model between the absolute value of the contrasts
#' and the ancestral state inferred at the node at which the contrast was calculated
#'
#' @param unit.tree a 'unit.tree' object
#'
#' @details This summary statistic is the estimate of the slope from fitting
#' a linear model between the contrasts and the inferred ancestral state at the node
#' at which the contrasts were calculated. It is used
#' to evaluate whether the model is adequately capturing variation relative to ancestral state.
#' The ancestral states are calculated using the contrast based approach and are not technically
#' ancestral state estimates (see Felsenstein 1985); this is conceptually and statistically
#' different from the maximum likelihood approach to reconstructing ancestral states.
#' The contrast based approach is
#' used here as there is a clear expectation of the relationship between the the contrasts and
#' trait value. If the generating model is correct, we expect there not to be a relationship
#' between the contrasts and their node height (i.e. slope ~ 0).
#'
#' The slope of the contrasts and their variances is included as a default summary statistic
#' in the function \code{\link{summ.stats}} 
#' 
#' @return the estimated slope parameter
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{def.summ.stats}}, \code{\link{stats::lm}}, \code{\link{ape::pic}}
#'
#'
#' @export slope.pic.asr
#'
#' @examples
#' data(geospiza)
#' td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
#' phy <- td$phy
#' dat <- td$data[,"wingL"]
#' unit.tree <- as.unit.tree(phy, dat)
#'
#' ## estimate slope
#' slope.pic.asr(unit.tree)
#' 
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










#' @title Internal function for getting default summary statistics
#'
#' @description Creates a list of the default summary statistics to be used
#' to asses model adequacy
#'
#' @details The following summary statistics are produced by this function:
#' \enumerate{
#'   \item{reml.sigsq: }{The mean of the squared contrasts. This is equivalent to the REML estimate of sigsq.}
#'
#'   \item{var.con: }{The variance of the absolute value of the contrasts.}
#'
#'   \item{slope.con.var: }{The slope of a linear model fit between the contrasts and their expected variances.}
#'
#'   \item{slope.con.asr: }{The slope of a linear model fit between the contrasts and their inferred ancestral state.}
#'
#'   \item{slope.con.nh: }{The slope of a linear model fit between the contrasts and the node height at which they were calculated.}
#'
#'   \item{ks.pic: }{The D-statistic from a KS test comparing the distribution of the contrasts to a normal distribution  with mean 0 and variance equal to the square root of the squared mean of the contrasts.}
#'  }
#'
#' @return a named list of functions
#'
#' @export def.summ.stats
#'
#' @keywords internal
#'
#' @seealso \code{\link{summ.stats}}, \code{\link{sigsq.reml}}, \code{\link{var.pic}}, \code{\link{slope.con.bl}}, \code{\link{slope.con.nh}}, \code{\link{slope.con.asr}}, \code{\link{ks.pic}}
#'
#'
#' @examples
#' ## get default summary stats
#' stats <- def.summ.stats
#' stats
def.summ.stats <- function()
    list("reml.sigsq"=sigsq.reml, "var.con"=var.pic, "slope.con.var"=slope.pic.var, "slope.con.asr"=slope.pic.asr, "slope.con.nh"=slope.pic.nh, "ks.dstat"=ks.pic)






#' @title Internal function to check summary stats
#'
#' @description Makes sure summary statistics sent to \code{\link{summ.stats}} are in the correct format
#'
#' @param stats a named list of summary statistics
#'
#' @return named list of summary statistics after being checked.
#'
#' @seealso \code{\link{summ.stats}}
#'
#' @export check.summ.stats
#'
#' @keywords internal
#'
#'
#' @examples
#' ## produce list of summary stats (using default summary statistics)
#' check.summ.stats(def.summ.stats())
#'
#' ## use custom list of summary statistics
#' foo <- function(x) mean(x$pics[,"contrasts"])
#'
#' my.stats <- list(mean.of.contrasts=foo)
#' check.summ.stats(my.stats)
#
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
