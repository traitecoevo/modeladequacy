## this is a script to run an analysis using a combination of old and new funcitonality

library(arbutus)
library(geiger)




## fxns from old code to parse fitContinuous result

parsefitContinuous <- function(fit){
	
	arg.names <- argn(fit$lik)[1]
	
	all.models <- c("BM", "OU", "EB", "lambda", "kappa", "delta", "trend")
	
	param <- c("sigsq", "alpha", "a", "lambda", "kappa", "delta", "slope")
	
	tmp <- cbind(all.models, param)
	
	model <- as.character(tmp[tmp[,"param"] == arg.names, "all.models"])
	
	sigsq <- fit$opt$sigsq
	
	z0 <- fit$opt$z0
	
	if (model == "BM"){ # BM does not have a model specific parameter
		model.p <- 0
	} else {
		model.p <- fit$opt[[1]]
	}
	
	if ("SE" %in% names(fit$opt)){
			SE <- fit$opt$SE
	} else {SE <- 0}

	par <- c(sigsq, z0, model.p, SE)
	names(par) <- c("sigsq", "z0", arg.names, "SE")
	
	return(list(model=model, par=par))
	
}




## function from old code  to rescale phylogeny from parsefitContinuous output

rescale.fitContinuous <- function(phy, fitCont){

    if (fitCont$model == "BM"){ ## if BM we don't need to transform the tree

		unit.tree <- phy
		unit.tree$edge.length <- unit.tree$edge.length * fitCont$par[1] 
		
		nt <- Ntip(unit.tree)
		tips <- unit.tree$edge[,2] <= nt
		unit.tree$edge.length[tips] <- unit.tree$edge.length[tips] + fitCont$par[4]^2
	}
	
	if (fitCont$model != "BM"){ ## for everything else
	
		unit.tree <- geiger:::rescale.phylo(phy, fitCont$model, fitCont$par[3])
		unit.tree$edge.length <- unit.tree$edge.length * fitCont$par[1]
		
		nt <- Ntip(unit.tree)
		tips <- unit.tree$edge[,2] <= nt
		unit.tree$edge.length[tips] <- unit.tree$edge.length[tips] + fitCont$par[4]^2
	
	}

    unit.tree
 }



## Now back to our regularly scheduled program...




## load in the geospiza data
data(geospiza)

tree <- geospiza$phy
data <- geospiza$dat[,1]


## fit BM  model with fitContinuous
fit.bm <- fitContinuous(phy=tree, dat=data, model="BM")



## USING OLD CODE
## get fitContinuous output
out  <- parsefitContinuous(fit.bm)

## USING OLD CODE
## rescale tree based on model parameters
rescaled.tree <- rescale.fitContinuous(tree, out)





## Using new code in arbutus
## Attach the data to the rescaled.tree; append data
unit.tree <- as.unit.tree(rescaled.tree, data)

## Calculate summary statistics on observed data set
ss.obs <- summStat(unit.tree)

## Simulate datasets on unit.tree
sims <- sim.charUnit(unit.tree=unit.tree, nsim=1000)

## Calculate summary statistics on simulated datasets
ss.sim <- summStat(sims)

## Compare observed to simulated summary statistics
res <- compare.summStat(summ.stats.obs = ss.obs, summ.stats.sim = ss.sim)

## look at results
res
