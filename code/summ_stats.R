## Model Adequacy fxns
## 
## Matthew Pennell
## 07/26/13

##
## libraries to require
library(geiger)
library(ggplot2)
library(multicore)


## Functions to make summary statistics
## Take a unit tree
## Return a fxn which takes tip data as an argument


## SS1: REML estimate of sigsq using phylogenetic independent contrasts
## Equal to the mean of the squared values

make.sigsqReml <- function(unit.tree){
	
	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class phylo"))
	}
	
	.sigsqReml <- function(data){
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data
	
	## Take pics
	pics <- pic(data, tree) 
	
	## Get mean of squared pics 
	remlss <- mean(pics^2)	
	
	remlss
	
	}
	
	list(fxn=.sigsqReml, unit.tree=unit.tree)
	
}





## SS2: Kolmorgorov-Smirnoff D statistic
## Compares distribution of pics with a normal distribution 
## (0, sqrt[mean contrasts^2])


make.ksPic <- function(unit.tree){

	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
	.ksPic <- function(data){
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data

	## Take pics
	pics <- pic(data, tree) 
	
	## SD is the mean of the squared contrasts
	sd <- sqrt(mean(pics^2))
	
	## Create null distribution 	
	nulldist <- rnorm(10000, mean=0, sd=sd)
	
	## KS test; return D statistic 	
	ksBM <- ks.test(pics, nulldist)$statistic
	 	
	as.numeric(ksBM)
	
	}
	
	list(fxn=.ksPic, unit.tree=unit.tree)
}





## SS3: Variance in pics
## Calculates variance in the absolute value of the pics

make.varPic <- function(unit.tree){

	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
	.varPic <- function(data){
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data

	## Take pics
	pics <- pic(data, tree) 
	
	## Get var of absolute value of pics
	var <- var(abs(pics)) 
	
	var

	}
	
	list(fxn=.varPic, unit.tree=unit.tree)
}






## SS4: Slope between absolute value of contrasts and their variances
## Tests for variance with respect to branch lengths

make.slopePicBl <- function(unit.tree){

	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
	.slopePicBl <- function(data){
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data

	## Take pics, get variance
	pics <- pic(data, tree, var.contrasts=TRUE)
	
	## Sqroot of the variance for each contrast
	sd.pic <- sqrt(pics[,"variance"]) 
	
	## Absolute value of the contrasts
	abs.pic <- abs(pics[,"contrasts"]) 
	
	## Fit linear model
	con.var <- lm(abs.pic ~ sd.pic) 
	
	## Extract slope
	slope <- con.var$coefficients["sd.pic"]
	
	if (is.na(slope) == TRUE){
		return(NA)
	}
	 
	as.numeric(slope)
	
	}
	
	list(fxn=.slopePicBl, unit.tree=unit.tree)

}





## SS5: Slope between inferred absolute value of contrasts and inferred ancestral state
## Test for variance with respect to state

make.slopePicAnc <- function(unit.tree){
	
	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
	.slopePicAnc <- function(data){
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data

	## Take absolute value of pics
	abs.pic <- abs(pic(data, tree))
	
	## Get ancestral states
	anc.st <- ace(data, tree, method="pic")$ace
	
	## Fit linear model
	con.anc <- lm(abs.pic ~ anc.st)
	
	## Extract slope
	slope <- con.anc$coefficients["anc.st"]
	
	if (is.na(slope) == TRUE){
		return(NA)
	}
	 
	as.numeric(slope)	
	
	}
	
	list(fxn=.slopePicAnc, unit.tree=unit.tree)
	
}




## SS6: Slope between absolute value of contrasts and the node height
## Test for variance with respect to node height
## Aka: Node-Height Test

make.slopePicNh <- function(unit.tree){
	
	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
	.slopePicNh <- function(data){
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data
	
	## take the absolute value of pics
	abs.pic <- abs(pic(data, tree))
	
	## get the node heights
	node.h <- branching.times(tree)
	
	con.nh <- lm(abs.pic ~ node.h)

	## Extract slope
	slope <- con.nh$coefficients["node.h"]
	
	if (is.na(slope) == TRUE){
		return(NA)
	}
	 
	as.numeric(slope)
	
	}
	
	list(fxn=.slopePicNh, unit.tree=unit.tree)
}








## Set default summary statistics

.defaultSS <- function(){
	dss <- c(make.sigsqReml, make.ksPic, make.varPic, make.slopePicBl, make.slopePicAnc, make.slopePicNh)
	
	dss.names <- c("REML.sigsq", "KS.Dstat", "Var.pic", "Slope.pic.var", "Slope.pic.anc", "Slope.pic.nh")
	
	list(f=dss, f.name=dss.names)
	
}






## Function to create default summary stat fxns

make.defStats <- function(unit.tree){

	dss <- .defaultSS()
	ss <- lapply(dss$f, function(f) f(unit.tree)) 
	names(ss) <- dss$f.name
	ss
	
}





## Calculate summary statistics for a given data set
## Takes two arguments: 
## 1. data -- a named vector of tip data
## 2. stats -- a list of fxns to compute (see make.defStats)


traitStat <- function(data, stats){

	ss <- lapply(stats, function(f) f$fxn(data))
	stats <- do.call(c, ss)
	stats
}






## Simulate data on unit tree
## Wrapper fxn for sim.char
## Takes unit.tree
## Takes nsim (specify number of sims to do)

sim.charUnit <- function(unit.tree, nsim){

	## Simulate trait data under BM with rate 1
	sims <- sim.char(phy=unit.tree, par=1, nsim=nsim, model="BM")
	
	## Convert format so that it can be used with lapply
	data <- lapply(c(1:nsim), function(x) sims[,,x])
	
	data

}





## Wrapper function which simulates n data sets and calculates summary statistics
## Trait Stat parametric bootstrap
## Takes the following as arguments
## Named list of functions (same as given to traitStat
## nsim
## Gets unit.tree from fxns

traitStatPb <- function(stats, nsim){

	## Collect unit trees
	ut <- lapply(stats, function(x) return(x$unit.tree))
	
	## Check to make sure that they are all identical
	if (length(unique(ut)) != 1){
		return(print("Not all statistic fxns were built with the same unit tree"))
		}

	## Simulate nsim data sets
	sim.data <- sim.charUnit(ut[[1]], nsim)
	
	## Calculate summary statistics across all datasets
	sim.ss <- lapply(sim.data, function(x) traitStat(x, stats))
	
	## Convert to data.frame
	sim.ss <- do.call(rbind, sim.ss)
	
	sim.ss

}





