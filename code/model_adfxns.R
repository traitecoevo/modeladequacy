## Model Adequacy Functions
##
## Written by Matt Pennell
## July 20, 2013
##
## Intro: this file contains all of the code required to evaluate model adequacy
## for continuous character evolution
##
##


library(geiger)
library(ggplot2)
library(phytools)





## Summary statistics

## These all take a unit.tree and the data as a named vector

## SS1: REML estimate of sigsq using phylogenetic independent contrasts
## Equal to the mean of the squared values

sigsqReml <- function(unit.tree, data){
	
	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class phylo"))
	}
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data
	
	## Take pics
	pics <- pic(data, tree) 
	
	## Get mean of squared pics 
	remlss <- mean(pics^2)	
	
	return(remlss)
	
}





## SS2: Kolmorgorov-Smirnoff D statistic
## Compares distribution of pics with a normal distribution 
## (0, sqrt[mean contrasts^2])


ksPic <- function(unit.tree, data){

	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
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
	 	
	return(ksBM)

}





## SS3: Variance in pics
## Calculates variance in the absolute value of the pics

varPic <- function(unit.tree, data){

	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data

	## Take pics
	pics <- pic(data, tree) 
	
	## Get var of absolute value of pics
	var <- var(abs(pics)) 
	
	return(var)

}






## SS4: Slope between absolute value of contrasts and their variances
## Tests for variance with respect to branch lengths

slopePicBl <- function(unit.tree, data){

	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
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
	 
	return(as.numeric(slope))	

}





## SS5: Slope between inferred absolute value of contrasts and inferred ancestral state
## Test for variance with respect to state

slopePicAnc <- function(unit.tree, data){
	
	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
	## check tree and data names
	td <- treedata(phy=unit.tree, data=data)
	tree <- td$phy
	data <- td$data

	## Take absolute value of pics
	abs.pic <- abs(pic(data, tree))
	
	## Get ancestral states
	anc.st <- fastAnc(tree, data)
	
	## Fit linear model
	con.anc <- lm(abs.pic ~ anc.st)
	
	## Extract slope
	slope <- con.anc$coefficients["anc.st"]
	
	if (is.na(slope) == TRUE){
		return(NA)
	}
	 
	return(as.numeric(slope))	
	
}




## SS6: Slope between absolute value of contrasts and the node height
## Test for variance with respect to node height
## Aka: Node-Height Test

slopePicNh <- function(unit.tree, data){
	
	## check to make tree is a "phylo" object
	if (class(unit.tree) != "phylo"){
		return(print("Unit.tree must be of class 'phylo'"))
	}
	
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
	 
	return(as.numeric(slope))	
	

}





## Wrapper function for calculating all SS1:SS6

traitSumm <- function(unit.tree, data){
	
	## create empty vector
	res <- vector()
	
	res[1] <- sigsqReml(unit.tree, data)
	res[2] <- ksPic(unit.tree, data)
	res[3] <- varPic(unit.tree, data)
	res[4] <- slopePicBl(unit.tree, data)
	res[5] <- slopePicAnc(unit.tree, data)
	res[6] <- slopePicNh(unit.tree, data)
	
	names(res) <- c("REML.sigsq", "KS.Dstat", "Var.pic", "Slope.pic.var", "Slope.pic.anc", "Slope.pic.nh")
	
	return(res)
	
}



