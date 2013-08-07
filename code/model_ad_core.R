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
library(multicore)




## Summary statistics

## These all take a unit.tree and the data as a named vector

## They return a function

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
	
	return(remlss)
	
	}
	
	return(.sigsqReml)
	
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
	anc.st <- ace(data, tree, method="pic")
	
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
## Outputs a named vector 6

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




## ML version of model adequacy

## Simulate n data sets on unit tree
## Under BM with a rate of 1
## Calculate summary statistics on each data set
## Built to utilize multiple cores with multicore package
## The ... represent options for the mclapply function used internally



## internal function for doing the simulation


.traitSummSim <- function(unit.tree){
		
	## simulate data
	data <- sim.char(unit.tree, par=1, model="BM")[,,]
		
	## get summary stats
	ss <- traitSumm(unit.tree, data)
		
	return(ss)
}





## Wrapper function which computes summary stats across simulated data sets


traitSummSim <- function(unit.tree, n.sims, n.cores=1, ...){
	
	res <- mclapply(c(1:n.sims), function(x) .traitSummSim, mc.cores=n.cores, ...)
	
	res <- do.call(rbind, res)
	return(res)
}




## Calculate p-values comparing observed to simulated data set (ML method)

traitSummPvalue <- function(summ.stats.data, summ.stats.sim){
	
	## check to make sure names are same
	if (names(summ.stats.data) != colnames(summ.stats.sim)){
		return(print("Summary stats from data to not match summary stats from simulation"))
	}
	
	p.values <- vector()
	
	## calculate p value for each summary statistics
	for (i in 1:length(summ.stats.data)){
		
		## empirical summary stat
		em <- summ.stats.data[i]
		
		## simulated summary stats
		sim <- as.vector(summ.stats.sim[,i])
		
		## calculate p-value
		## one-sided test
		p.1 <- length(which(sim >= em))/(length(sim) + 1)
		p.2 <- length(which(sim < em))/(length(sim) + 1)
		p <- min(p.1, p.2)
		
		## bind to other stats
		p.values <- c(p.value, p)
	}
	
	names(p.values) <- names(summ.stats.data)
	return(p.values)
}



