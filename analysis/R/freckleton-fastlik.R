require( MCMCpack)
require(mcmc)


# This code is intended to show some side-by-side examples of algorithms for fitting 
# phylogenetic models with contrasts, compared with the method of direct matrix
# inversion. It is not a fully fledged R-package: hopefully what it will serve to do
# is to act as a benchmark  / example for those who have already developed R 
# packages for fitting comparative models. 



# ----------------------------------------------------------------------
# Compute a likelihood for given parametets & data
# ----------------------------------------------------------------------
likGeneral  <- function(x, V, mu, s2, lambda = 1) {
	
	# Calculates the likelihood
	lik <- function(x, V, mu, s2) {
		n <- dim(x)[1]
		k <- dim(x)[2] 
		xT <- c(t(x))
		ks2V <- kronecker( V, s2)
		ks2Vi <- solve( ks2V)
		X <- rep( mu, dim(x)[1])
		res <- crossprod( ( xT - X), ks2Vi %*% (xT - X) ) 
		ll <- -0.5 * ( n * k * log( 2 * pi ) + determinant( ks2V, log = TRUE)$modulus[1]  + res )
		return( list(ll = ll, mu = mu, s2 = s2) )
		}
			
	x <- orderD(x)
	if(is.vector(x)  == TRUE) x <- matrix(x, ncol  =1) 
	V <- orderV(V)
	Vmat <- lamTrans(V, lambda)
	return(lik(x, Vmat, mu, s2))
	}

# ----------------------------------------------------------------------
# Lambda transformation - applied to vcv matrix
# ----------------------------------------------------------------------
lamTrans <- function(V, lambda) {
	V1 <- V
	diag(V1) <- 0
	V2 <- diag( diag(V), ncol = length(V[1,]), nrow = length(V[,1]))
	Vmat <- V1 * lambda + V2
	return(Vmat)
	}

# ----------------------------------------------------------------------
# Compute phylogenetically corrected mean vector
# ----------------------------------------------------------------------
estMean <- function(y, V) {
	iV <- solve(V, tol = .Machine$double.eps)
	xdum <- matrix(rep(1, dim(y)[1] ) , nrow = dim(y)[1])
	xVix <- crossprod(xdum, iV %*% xdum)
	xViy <- crossprod(xdum, iV %*% y)
	mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy # This is a bad thing!
	return(mu)
	}
	
# ----------------------------------------------------------------------
# Compute phylogenetically corrected variance matrix for traits
# ----------------------------------------------------------------------
estVar <- function(x, V) {
	mu <- estMean(x, V)
	iV <- solve(V, tol = .Machine$double.eps)
	e <- x
	for( i in 1:length(mu) ) e[,i] <- e[,i]  - mu[i]
	s2 <- crossprod(e, iV%*%e)
	n <- dim(x)[1]
	k <- dim(x)[2] 
	return(s2 / n)
	}

# ----------------------------------------------------------------------
# This section of code fits lambda to data on multiple traits
# Using maximum likelihood via direct solution of the likelihood
# ----------------------------------------------------------------------
likLambda <- function(x, V, lambda) {

	# Computes the likelihood	 
	lik <- function(x, V) {
		mu <- estMean(x, V)
		s2 <- estVar(x, V)
		n <- dim(x)[1]
		k <- dim(x)[2] 
		xT <- c(t(x))
		ks2V <- kronecker( V, s2)
		ks2Vi <- solve( ks2V)		
		X <- rep( mu, dim(x)[1])
		res <- crossprod( ( xT - X), ks2Vi %*% (xT - X) ) 
		ll <- -0.5 * ( n * k * log( 2 * pi ) + determinant( ks2V, log = TRUE)$modulus + res )
		return( list(ll = ll, mu = mu, s2 = s2) )
		}
		
	x <- as.matrix(x)
	x <- orderD(x)
	V <- orderV(V)
	if(is.vector(x)  == TRUE) x <- matrix(x, ncol  =1) 
	Vmat <- lamTrans(V, lambda)
	return(lik(x, Vmat))
	}

# -----------------------------------------------
# Maximizes likelihood using direct method
# -----------------------------------------------
maxLikLambda <- function(x, V) {
		ll <- function(lambda) return(-1* likLambda(x, V, lambda)$ll )
		optimize(ll, c(0,1)) 
		}

# ----------------------------------------------------------------------
# This section of code fits multiple lambdas to data on multiple traits
# Using maximum likelihood via direct solution of the likelihood
# ----------------------------------------------------------------------
likMultilambda <- function(x, phylomat, lambdas) {

	lik <- function(x, V) {

		n <- dim(x)[1]
		k <- dim(x)[2] 
		xT <- c(t(x))
		
				
		mu <- rep(0, k)

		for( i in 1:k) {
			Vmat <- lamTrans(V, lambdas[i]) 
			mu[i] <- estMean( matrix( x[,i], n, 1), Vmat)
			}
		
		s2 <- matrix(0, k, k)
		
		for( i in 1:k) {
			for( j in 1:k) {
					Vmati <- lamTrans(V, lambdas[i] )
					Vmatj <- lamTrans(V, lambdas[j] )
					#Vmat <- sqrt( Vmati) * sqrt(Vmatj)	
					Vmat <- lamTrans(V, sqrt( lambdas[i] * lambdas[j]) )
				
					mui <- estMean( matrix( x[,i], n, 1), Vmat)
					muj <- estMean( matrix( x[,i], n, 1), Vmat)
					
					xi <- matrix( x[,i] - mui, n, 1)
					xj <- matrix( x[,j] - muj, n, 1)
					s2[i,j] <- crossprod( xj, solve(Vmat) %*% xi)
					}
				}
		
		s2 <- s2 / n
		
		ks2V <- matrix(0, nrow  = k*n, ncol = k*n)

			for(j in 1:n) {
				for(m in 1:n) {
					for(i in 1:k) {
					for(l in 1:k) {
					r <- i + (j - 1) * k 
					c <- l + (m - 1) * k
					if( j != m ) { ks2V[r,c] <- s2[i, l] * sqrt(lambdas[i] * lambdas[l]) * V[j, m] } else{ ks2V[r,c] <- s2[i, l] * V[j, m] }
					}
				}
			}
		}

		
		invVmat <- solve( ks2V )
		
		X <- rep( mu, dim(x)[1])
		res <- crossprod( ( xT - X), invVmat %*% (xT - X) ) 
		ll <- -0.5 * ( n * k * log( 2 * pi ) + log( det( ks2V) ) + res )
		
		return( list(ll = ll, mu = mu, s2 = s2) )
		}
	
	x <- as.matrix(x)
	x <- orderD(x)
	V <- orderV(phylomat)

	return(lik(x, V))
	}



# ----------------------------------------------------------------------
# This maxmimzes the likelihood using the direct method
# ----------------------------------------------------------------------
maxLikMultilambdas <- function(x, V) {
	ll <- function(lambdas) return( -1 * likMultilambda(x, V, lambdas)$ll )
	op <- optim( rep(0.5, dim(x)[2]), ll, method = "L-BFGS-B", lower = rep(0, dim(x)[2]), upper = rep(1, dim(x)[2]))
	
	lambdas <- op$par
	names(lambdas) <- colnames(x)
	ML <- op$value
	MLmod <- likMultilambda( x, V, lambdas) 
	
	return( list( ML = ML, lambda = lambdas, MLmod = MLmod) )

	}


# ----------------------------------------------------------------------
# Samples from the posterior using MCMC sampling
# ----------------------------------------------------------------------
mcmcLik <- function(x, V, scale = 1, nbatch = nbatch, burnin = 1000, priorMu = 0, priorVar = 0) {
	
	k <- dim(x)[2]
	nvars <- k^2 - (k^2 - k ) / 2

	likFunc <- function(pars) {
		mu <- pars[1:k]
		s2 <- matrix(0, nrow = k, ncol = k)
		idx <- k+1
		for( i in 1:k) {
			for( j in i:k) {
				s2[i, j] <- pars[idx]
				s2[j, i] <- pars[idx]
				idx <- idx +1
				}
			}
		for( i in 1:k) {
			for( j in i:k) {
				if( i != j) {
					
					s <- exp(  s2[i,j] ) / ( 1 + exp( s2[i,j] ) )
					s <- ( s * 2 ) - 1
					
					s2[i, j] <-  s * sqrt(s2[i,i] * s2[j,j])
					s2[j, i] <-  s * sqrt(s2[i,i] * s2[j,j])
				 	}
				}
			}
		
		
		priorMu <- sum( dnorm(mu, 0, 4, log = TRUE) )
	#	priorS2 <- log( diwish(s2, 4, diag(1, k,k) ) )
		b <- k + 2
		z <- b - k -1
		
		priorSigma <-   z * diag(1, k) 
		
		priorS2 <- log( diwish(s2, b,  priorSigma) )		
		L <- likGeneral(x, V, mu, s2)$ll  + priorMu + priorS2
		
		
		return( L )
		}
	
	ml <- likLambda(x, V, 1)
	s2 <- ml$s2
	idx <- 1
	inits <- rep(0, nvars - k)
	idx <- 1
	for( i in 1:k) {
		for( j in i:k) {
			if( i != j) {					
				s2[i,j] <- s2[i,j] / sqrt( s2[i,i] * s2[j,j] )
				s2[j,i] <- s2[j,i] / sqrt( s2[i,i] * s2[j,j] )
				s2[i,j ]<- ( s2[i,j] + 1 ) / 2
				inits[idx] <- log( s2[i,j] /  (1 - s2[i,j] ))
			} else { inits[idx] <- s2[i,j]}
			
			idx <- idx +1
		}
	}
	
	init <- c(ml$mu, inits)
	burnin <- metrop( likFunc, initial = init, nbatch = burnin, scale = scale)
	metrop( burnin, nbatch =nbatch)

	}



# ---------------------------------------------------------------------------------
# --------- Methods based on using contrasts for fast likelihood calculations
# ---------------------------------------------------------------------------------

# -----------------------------------------------
# ------------- Lambda Transformation ------
# -----------------------------------------------
lambdaTrans <- function( phy, lambda, height = NULL ) {
    n <- length(phy$tip.label)
    if( is.null( height) == TRUE) height <- tree.height(phy)
    interns <- which( phy$edge[, 2] > n )
    externs <- which(  phy$edge[, 2] <= n)
    phy$edge.length[ interns ] <-  phy$edge.length[ interns ]  * lambda
    phy$edge.length[ externs ] <-  phy$edge.length[ externs ]  * lambda + (1 - lambda) * height
	return( phy )
	}

# -----------------------------------------------
# ------------- Stolen from APE --------------
# -----------------------------------------------
tree.height <- function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\".")
    if (is.null(phy$edge.length)) 
        stop("the tree has no branch lengths.")
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in 1:dim(phy$edge)[1]) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 
        1]] + phy$edge.length[i]
 	return(xx[1:n])
}


# -----------------------------------------------
# ----- My modification to Emmanuel's code 
# -----------------------------------------------
picRob <- function (x, phy, scaled = TRUE, var.contrasts = TRUE) 
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



clikGeneral <- function(x, phylo, mu, s2, lambda = 1) {
	if(lambda != 1) phylo <- lambdaTrans(phylo, lambda)
	n <- dim(x)[1]
	k <- dim(x)[2]
	picx <- apply( x, 2, function(i) { names(i) <- rownames(x); picRob(i, phylo)} )
	cx <-lapply(picx,  function(x) x$contr[,1] )
	cx <- matrix( unlist(cx), ncol = k)
	W <- as.matrix(s2)
	iW <- solve(W)
	vars <- c( picx[[1]]$contr[,2], picx[[1]]$V )
	
	muRoot <- lapply( picx, function(x) x$root.state)
	muRoot <- unlist(muRoot)
	cRoot <- ( muRoot - mu ) / sqrt( picx[[1]]$V )
	S <- 0
	for(i in 1:(n-1) ) S <- S + cx[i,] %*% iW %*% cx[i,]
	S <- S + cRoot %*% iW %*% cRoot
	LL <- -0.5 * ( n * k * log( 2 * pi) + k * sum(log(vars) ) + n * log( det(W ) ) + S )
	
	return( list( LL = LL, mu = mu, s2 = W) )
	}

# ----------------------------------------------------------------------
# This section of code fits the model
# Using maximum likelihood via the fast contrast method
# ----------------------------------------------------------------------
clikLambda <- function( x, phylo, lambda) {
	if(lambda != 1) phylo <- lambdaTrans(phylo, lambda)
	n <- dim(x)[1]
	k <- dim(x)[2]
	picx <- apply( x, 2, function(i) { names(i) <- rownames(x); picRob(i, phylo)} )
	cx <-lapply(picx,  function(x) x$contr[,1] )
	cx <- matrix( unlist(cx), ncol = k)
	W <- crossprod(cx, cx) / n
	iW <- solve(W)
	vars <- c( picx[[1]]$contr[,2], picx[[1]]$V )
	
	S <- 0
	for(i in 1:(n-1) ) S <- S + cx[i,] %*% iW %*% cx[i,]
	LL <- -0.5 * ( n * k * log( 2 * pi) + k * sum(log(vars) ) + n * log( det(W ) ) + S )
	
	mu <- lapply( picx, function(x) x$root.state)
	mu <- unlist(mu)
	return( list( LL = LL, mu = mu, s2 = W) )
	}
	
# ----------------------------------------------------------------------
# This section of code fits the model
# Using REML via the fast contrast method
# ----------------------------------------------------------------------
REMLclikLambda<- function( x, phylo, lambda) {
	if(lambda != 1) phylo <- lambdaTrans(phylo, lambda)
	n <- dim(x)[1]
	k <- dim(x)[2]
	picx <- apply( x, 2, function(i) { names(i) <- rownames(x); picRob(i, phylo)} )
	cx <-lapply(picx,  function(x) x$contr[,1] )
	cx <- matrix( unlist(cx), ncol = k)
	W <- crossprod(cx, cx) / ( n - k ) 
	iW <- solve(W)
	vars <- c( picx[[1]]$contr[,2]) #, picx[[1]]$V )
	
	S <- 0
	for(i in 1:(n-1) ) S <- S + cx[i,] %*% iW %*% cx[i,]
	LL <- -0.5 * ( (n - k ) * k * log( 2 * pi) + k * sum(log(vars) ) + (n - k) * log( det(W ) ) + S )
	
	mu <- lapply( picx, function(x) x$root.state)
	mu <- unlist(mu)
	return( list( LL = LL, mu = mu, s2 = W) )
	}

# ----------------------------------------------------------------------
# This maxmimzes the likelihood using the contrast methods
# ----------------------------------------------------------------------
maxClikLambda <- function(x, phylo) {
		ll <- function(lambda) return(- 1* clikLambda(x, phylo, lambda)$LL )
		op <- optimize(ll, c(0 , 1))
		lambda <- op$minimum
		ML <- op$objective
		MLmod <- clikLambda( x, phylo, lambda)
		return( list( ML = ML, lambda = lambda, MLmod = MLmod) )	 
	}
	
# ----------------------------------------------------------------------
# This section of code fits multiple lambdas to data on multiple traits
# Using maximum likelihood via the fast contrast method
# ----------------------------------------------------------------------
clikMultilambda <- function( x, phylo, lambdas) {
	n <- dim(x)[1]
	k <- dim(x)[2]
	
	picx <- list()
	
	for(i in 1:k) {
		xi <- x[,i]
		names(xi) <-  rownames(x)
		picx[[i]] <- picRob( xi, lambdaTrans(phylo, lambdas[i] ) )
		}
	
	cx <-lapply(picx,  function(x) x$contr[,1] )
	cx <- matrix( unlist(cx), ncol = k)
	W <- crossprod(cx, cx) / n
	iW <- solve(W)
	
	vars <- NULL 
	for( i in 1:k) vars <- c(vars, picx[[i]]$contr[,2], picx[[i]]$V ) 
	
	S <- 0
	for(i in 1:(n-1) ) S <- S + cx[i,] %*% iW %*% cx[i,]
	LL <- -0.5 * ( n * k * log( 2 * pi) +  sum(log(vars) ) + n * log( det( W ) ) + S )
	mu <- lapply( picx, function(x) x$root.state)
	mu <- unlist(mu)
	return( list( LL = LL, mu = mu, s2 = W) )
	}


# ----------------------------------------------------------------------
# This section of code fits multiple lambdas to data on multiple traits
# Using maximum likelihood via the fast contrast method
# ----------------------------------------------------------------------
clik2Multilambda <- function( x, phylo, lambdas) {
	n <- dim(x)[1]
	k <- dim(x)[2]
	picx <- list()
	
	for(i in 1:k) {
		xi <- x[,i]
		names(xi) <-  rownames(x)
		picx[[i]] <- pic.Rob( xi, lambdaTrans(phylo, lambdas[i] ) )
		}
		
	cx <-lapply(picx,  function(x) x$contr[,1] )
	cx <- matrix( unlist(cx), ncol = k)
	W <- matrix(0, k, k)
	
	for(j in 1:k){
 		for(i in 1:k) {
			xi <- x[,i]
			xj <- x[,j]
			names(xi) <-  rownames(x)
			names(xj) <-  rownames(x)
			lambda <- ( sqrt(lambdas[i]) * sqrt(lambdas[j]) )
			picxi<- pic.Rob( xi, lambdaTrans(phylo, lambda ) )$contr[,1]
			picxj<- pic.Rob( xj, lambdaTrans(phylo, lambda  ) )$contr[,1]
			W[i, j] <- crossprod( picxi, picxj)  / (n - k)
			}
		}

	iW <- solve(W)
	vars <- NULL 
	for( i in 1:k) vars <- c(vars, picx[[i]]$contr[,2], picx[[i]]$V ) 
	S <- 0
	for(i in 1:(n-1) ) S <- S + cx[i,] %*% iW %*% cx[i,]
	LL <- -0.5 * ( n * k * log( 2 * pi) +  sum(log(vars) ) + n * log( det( W ) ) + S )
	
	mu <- lapply( picx, function(x) x$root.state)
	mu <- unlist(mu)
	return( list( LL = LL, mu = mu, s2 = W) )
	}

# ----------------------------------------------------------------------
# This maxmimzes the likelihood using the contrast methods
# ----------------------------------------------------------------------
maxClikMultilambdas <- function(x, phylo) {
	ll <- function(lambdas) return( -1 * clikMultilambda(x, phylo, lambdas)$LL )
	op <- optim( rep(0.5, dim(x)[2]), ll, method = "L-BFGS-B", lower = rep(0, dim(x)[2]), upper = rep(1, dim(x)[2]) )
	lambdas <- op$par
	names(lambdas) <- colnames(x)
	ML <- op$value
	MLmod <- clikMultilambda( x, phylo, lambdas) 
	return( list( ML = ML, lambda = lambdas, MLmod = MLmod) )

	}

# ----------------------------------------------------------------------
# Samples from the posterior using MCMC sampling - uses contrasts
# ----------------------------------------------------------------------
mcmcLikContrasts <- function(x, phy, scale = 1, nbatch = 1000, burnin = 1000, priorMu = 0, priorVar = 0) {
	
	k <- dim(x)[2]
	nvars <- k^2 - (k^2 - k ) / 2

	likFunc <- function(pars) {
		mu <- pars[1:k]
		s2 <- matrix(0, nrow = k, ncol = k)
		idx <- k+1
		for( i in 1:k) {
			for( j in i:k) {
				s2[i, j] <- pars[idx]
				s2[j, i] <- pars[idx]
				idx <- idx +1
				}
			}
	
		
		priorMu <- sum( dnorm(mu, 0, 4, log = TRUE) )
#		priorS2 <- log( diwish(s2, 4, diag(1, k,k) ) )
		b <- k + 2
		z <- b - k -1
		
		priorSigma <-   z * diag(1, k) 
		
		priorS2 <- log( diwish(s2, b,  priorSigma) )
		
		L <- -1000
		
	
		if( det(s2) > 0) L <- clikGeneral(x, phy, mu, s2)$LL  + priorMu + priorS2
		
		
		if(abs(L) == Inf) L <- -1000
		
		return( L )
		}
	
	ml <- clikLambda(x, phy, 1)
	s2 <- ml$s2
	idx <- 1
	inits <- rep(0, nvars - k)
	idx <- 1
	for( i in 1:k) {
		for( j in i:k) {
			if( i != j) {					
				inits[idx] <- s2[i,j] 
				} else { inits[idx] <- s2[i,j]}
			
			idx <- idx +1
		}
	}
	
	init <- c(ml$mu, inits)
	burnin <- metrop( likFunc, initial = init, nbatch = burnin, scale = scale, debug = TRUE)
	metrop( burnin, nbatch =nbatch)
	}
	

# ----------------------------------------------------------------------
# fits PGLM models using contrasts
# ----------------------------------------------------------------------
	
	
# -----------------------------------------------
# --------PGLM based on contrasts -----------
# -----------------------------------------------
cpglm <- function( formula, data, phylo) {
	# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
		phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	
	cpic <- function(x) pic(x, phylo)
	cdata <- apply( data, 2, cpic)
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata, y = TRUE)
	
	# Calculate the log-likelihood
	vars <- pic( data[,1], phylo, var = TRUE)[,2]
	V <- picRob(data[,1], phylo)$V
	u <- residuals(fm) 
	n <- length(u)
	sigma2 <- sum(u^2) / (n + 1)
	logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars )) + sum(u^2) / sigma2  )
	logLik <- logLik - 0.5 * log(V)
	
	return(list( model = fm, logLik = logLik) )
	
	}
	


# -----------------------------------------------
# -- PGLM based on contrasts with lambda -
# -----------------------------------------------	
cpglmLambda <- function( formula, data, phylo, lambda) {

# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
		phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	
	phy <- lambdaTrans( phylo, lambda )
	cpic <- function(x) pic(x, phy)
	cdata <- apply( data, 2, cpic)
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata, y = TRUE)
	
	# Calculate the log-likelihood
	vars <- pic( data[,1], phy, var = TRUE)[,2]
	V <- picRob(data[,1], phy)$V	
	u <- residuals(fm) 
	n <- length(u)
	sigma2 <- sum(u^2) / n
	logLik <- -0.5 *( (n + 1) * log( 2 * pi * sigma2) + sum(log( vars )) + sum(u^2) / sigma2  )
	logLik <- logLik - 0.5 * log(V)
	
	return(list( model = fm, logLik = logLik) )
	
	}


	
# -----------------------------------------------
# -- PGLM based on contrasts with lambda -
# -- With REML -----------------------------
# -----------------------------------------------	
REMLcpglmLambda <- function( formula, data, phylo, lambda) {

# Prune data & tree
	idx.retain <- which( complete.cases( data ) == TRUE)
	idx.drop <- which( complete.cases( data ) == FALSE)
	if(length(idx.drop) >0){
		phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
	data <- data[idx.retain,]
	
	phy <- lambda.trans( phylo, lambda )
	cpic <- function(x) pic(x, phy)
	cdata <- apply( data, 2, cpic)
	cdata <- data.frame(cdata)
	colnames(cdata) <- colnames(data)
	fm <- lm(formula, cdata, y = TRUE, x = TRUE)
	k <- length( coef(fm) ) 
	# Calculate the log-likelihood
	vars <- pic( data[,1], phy, var = TRUE)[,2] # do not include root variance for REML
	V <- picRob(data[,1], phy)$V	
	u <- residuals(fm) 
	n <- length(u) 
	sigma2 <- sum(u^2) / (n - k )
	#sigma2 <- var(u) 
	
	logLik.check <- logLik( fm, REML = TRUE) - 0.5 *  sum(log( vars )) # works for checking
	
	
	logLik <- -0.5 *(  (n - k) * log( 2 * pi * sigma2) + sum(log( vars )) + (n-k)  + sum(log( det( crossprod(fm$x, fm$x) )) ))
	
	
	return(list( model = fm, logLik = logLik, logLik.check = logLik.check, s2 = sigma2, n = n, k = k, ssq= sum(u^2),  v = sum(log( vars ))) )
	
	}
	
	
# -----------------------------------------------
# -- PGLM based on contrasts ---------------
# ---- Estimation of lambda ------------------
cpglmEstlambda <- function( formula, data, phylo ) {
	
			optimizeLambda <- function( formula, data, phylo ) {
				optFun <- function(pars) -1 * cpglmLambda( formula, data, phylo, lambda = pars[1])$logLik
				of <- optim( c(0.5),  method = "L-BFGS-B", optFun, lower = c(0), upper = c(1) )
				return(of)
				}
		
		# Prune data & tree
		idx.retain <- which( complete.cases( data ) == TRUE)
		idx.drop <- which( complete.cases( data ) == FALSE)
		if(length(idx.drop) >0){
			phylo <- drop.tip( phylo, rownames(data)[ idx.drop ])}
		data <- data[idx.retain,]
		
		opt <- optimizeLambda( formula, data, phylo )
		MLlambda <- opt$par[1]
		ML <-  -1 * opt$value
		
		fm <- cpglmLambda( formula, data, phylo, MLlambda)$model
		
		return( list(model = fm, lambda = MLlambda, ML = ML)  )
	
	}

	
	
	
	
		
# ----------------------------------------------------------------------
# Functions for getting the phylo matrix and data in alphabetical order
#
orderV <- function(Vmat) {
	nms <- rownames(Vmat)
	snms <- sort(nms, index.return = TRUE)
	Vmat <- Vmat[snms$ix, snms$ix]

	return(Vmat)
	}

orderD <- function(dat) {	
	nms <- rownames(dat)
	snms <- sort(nms, index.return = TRUE)
	dat <- dat[snms$ix,]
	
	return(dat)
	}

sortVec <- function(data) {
	nms <- names(data)
	snms <- sort( nms, index.return = TRUE)
	sdata <- data[snms$ix] 
	names(sdata)  <- snms$x
	return(sdata)
}

