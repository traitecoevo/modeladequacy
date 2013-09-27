## these are the functions for calculating summary statistics
##
## Written by Matt Pennell
## Last modified 18/09/13



## arbutus:::summStat

## Fxn for calculating summary statistics across a unit.tree or a set of unit.trees

## 2 arguments
## unit.tree -- can be either a single unit.tree or a list of unit.trees

## summ.stats -- a named list of summary statistics. Default is NULL, meaning that the fxn uses the set of 6 default summary statistics


## returns a named data.frame with all of the summary statistics

## Note: data from unit.trees will be used

summStat <- function(unit.tree, stats=NULL){
  stats <- checkSummStats(stats)
  if (inherits(unit.tree, "unit.tree")){
    # One unit.tree
    res <- do.call(cbind.data.frame,
                   lapply(stats, function(x) x(unit.tree)))
  } else if (is.list(unit.tree) && length(unit.tree) > 0) {
    # Multiple unit.trees
    res <- do.call(rbind, lapply(unit.tree, summStat, stats))
  } else {
    # Anything else
    stop("'unit.tree' must be a single unit.tree, or a list of them")
  }
  res
}

## arbutus:::sigsqReml

## A function for computing the REML estimate of sigma^2
## Equal to the mean of the squared contrasts
## Takes a unit.tree object
## If data == NULL, uses the data supplied by unit.tree
## If data != NULL, creates unit.tree with input data

sigsqReml <- function(unit.tree, data=NULL){
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	mean(unit.tree$pics[,"contrasts"]^2)
}





## arbutus:::ksPic

## A function for comparing the distribution of the contrasts
## to a normal distribution with mean 0 and sd = sqroot(mean(pics^2))
## uses D-statistic from KS test

ksPic <- function(unit.tree, data=NULL){
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	## simulate a null distribution
	sd <- sqrt(mean(unit.tree$pics[,"contrasts"]^2))
	nd <- rnorm(10000, mean=0, sd=sd)

	## KS test
	ksbm <- ks.test(unit.tree$pics[,"contrasts"], nd)$statistic

	unname(ksbm)
}







## arbutus:::varPic

## A fxn for computing the variance in pic estimates
## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data

varPic <- function(unit.tree, data=NULL) {
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

        var(abs(unit.tree$pics[,"contrasts"]))
}








## arbutus:::slopePicBl

## Fxn for computing the slope between the contrasts and their standard deviations
## Used as a measure of variation with respect to branch length

## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data

slopePicBl <- function(unit.tree, data=NULL){
	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	## sd of the pics
	sd.pic <- sqrt(unit.tree$pics[,"variance"])

	## absolute value of the pics
	abs.pic <- abs(unit.tree$pics[,"contrasts"])

	## fit a linear model
	cv <- lm(abs.pic ~ sd.pic)

	unname(coef(cv)["sd.pic"])
}









## arbutus:::slopePicNh

## Fxn for computing the slope between the contrasts and the height of the node at which the contrast is computed
## Used as a measure of variation with respect to time

## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data

# NOTES:
##
## 1. Need to swap out branching.times() with a function that works for non-ultrametric trees (working on the fxn nodeAge to do this)

slopePicNh <- function(unit.tree, data=NULL){

	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	## node heights (see note 1)
	nh <- branching.times(unit.tree$phy)

	## absolute value of the pics
	abs.pic <- abs(unit.tree$pics[,"contrasts"])

	## fit a linear model

	ch <- lm(abs.pic ~ nh)

	unname(coef(ch)["nh"])
}








## arbutus:::slopePicAsr

## Fxn for computing the slope between the contrasts and the inferred ancestral state at the node
## Used as a measure of variation with respect to trait value

## Takes a unit.tree object
## If data == NULL, uses the data attached to unit.tree
## If data != NULL, creates unit.tree with input data



## NOTES:
##
## 1. Use "pic" method in ape. Not sure if this is appropriate or not.


slopePicAsr <- function(unit.tree, data=NULL){

	## make sure the unit.tree is of class 'unit.tree'
	## if not create a unit.tree using supplied data
	if (!inherits(unit.tree, "unit.tree"))
		unit.tree <- as.unit.tree(unit.tree, data)

	## reconstruct ancestors using pic method
	asr <- ace(unit.tree$data, unit.tree$phy, method="pic")$ace

	## absolute value of pics
	abs.pic <- abs(unit.tree$pics[,"contrasts"])

	## fit linear model
	ca <- lm(abs.pic ~ asr)

	unname(coef(ca)["asr"])
}









## arbutus:::defSummStats

## function to build list of default summary statistics
## used internally in traitStats

## takes no arguments

defSummStats <- function()
    list("reml.sigsq"=sigsqReml, "var.con"=varPic, "slope.con.var"=slopePicBl, "slope.con.asr"=slopePicAsr, "slope.con.nh"=slopePicNh, "ks.dstat"=ksPic)

checkSummStats <- function(stats) {
  if (is.null(stats))
    stats <- defSummStats()
  else if (!is.list(stats) || is.null(names(stats)) ||
           any(names(stats) == "" || length(stats) == 0))
    ## (note that the second check does not always enforce
    ## length(stats) > 0)
    stop("'stats' must be a named list with at least one element")
  else if (!all(sapply(stats, is.function)))
    stop("All elements of stats must be functions")
  stats
}
