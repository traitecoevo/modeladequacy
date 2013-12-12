library(geiger)
library(diversitree)


## read in data
phy <- read.tree("salicaceae.tre")
d <- read.csv("salicaceae.csv", as.is=TRUE)
states <- d[,"V1"]
names(states) <- d[,"X"]

SE <- 0.1039405



## fit OU using fitContinuous
fitc <- fitContinuous(phy, states, SE=SE, model="OU")

fitc$opt[c("sigsq", "alpha", "z0")]


## my result

#$sigsq
#[1] 3.041565e-17

#$alpha
#[1] 1.741193

#$z0
#[1] 2.064679


## fit OU using diveristree
lik <- make.ou(phy, states, states.sd=SE,
               control=list(method="pruning", backend="C"))
dt <- find.mle(lik, x.init=c(1,0.1, mean(states)))

dt$par

## my result


#           s2         alpha         theta
# 1.155839e-11  2.273490e+00 -8.623042e+03

lik(unlist(fitc$opt[c("sigsq", "alpha", "z0")]))
lik(coef(dt))

## The models are fitting fundamentally the same thing:
fitc$lik(coef(fitc)[1:2])
lik(unlist(fitc$opt[c("sigsq", "alpha", "z0")]))

## Profile likelihood over twice the observed range:
f <- function(x) {
  constraint <- as.formula(paste("theta ~", x))
  lik.theta <- constrain(lik, constraint)
  fit <- find.mle(lik.theta, x.init=c(1, .1))
  c(p=fit$lnLik, coef(fit))
}
pp <- seq(coef(fitc)[["z0"]], 2*coef(dt)[["theta"]], length.out=101)
ans <- sapply(pp, f)
plot(pp, ans["p",], type="o", pch=19, cex=.5, xlab="Theta",
     ylab="Log likelihood")

## So it increases in an unbounded way towards the left.  Weird.

## OK, so looking at:
## ```coffee
## attr(lik(coef(dt), intermediates=TRUE), "vals")
## ```
## gies the mean, variance and normalisation constant for the gaussian
## function at the root of the tree.  The mean is 2.621384e+81, which
## falls a little outside the range of observed values!

## But then so does the same for geiger:
## ```
## p.geiger <- unlist(fitc$opt[c("sigsq", "alpha", "z0")])
## attr(lik(p.geiger, intermediates=TRUE), "vals")
## ```
## which has a mean of 5.6e48.

## Notice the variances on these though; both right around
## 1.246727e+109.  That's the weird thing about OU; when you look
## *backwards* at OU (that is the distribution of where you came from,
## conditional on where you are), it looks completelty different to
## the forward time realisation of the process; the difference between
## the PDF mean and the OU mean grows exponentially with time, so that
## any perturbation from the OU mean is amplified.  But the variance
## grows exponentially too; this just reflects the memory erasing
## property of OU - you really could have come from basically
## anywhere, but it makes a little more sense that you have not
## actually crossed over the optimum.  This is the reason why there is
## no QuaSSE+OU model, because modelling the required state space is
## really really hard.

## So, back to the original problem.  This is using the VCV approach.
lik.vcv <- make.ou(phy, states, states.sd=SE,
                   control=list(method="vcv"))

## In contrast with the pruning approach, we can't see any difference
## with changing theta here, but that's only because we ignore it.
lik(coef(dt))     # 38.538
lik.vcv(coef(dt)) # 37.987

lik(p.geiger)     # 37.987
lik.vcv(p.geiger) # 37.987

lik2 <- make.ou(phy, states,
               control=list(method="pruning", backend="C"))
dt2 <- find.mle(lik2, x.init=c(1,0.1, mean(states)))

fitc2 <- fitContinuous(phy, states, model="OU")

## These look much closer:
dt2$lnLik     # 39.88955
fitc2$opt$lnL # 39.8895

## And these are closer but still not the same around the optimum
coef(dt2)
##        s2     alpha     theta
## 0.0194262 1.3878969 1.5721008

unlist(fitc2$opt[c("sigsq", "alpha", "z0")])
##      sigsq      alpha         z0
## 0.01942321 1.38774506 2.06620538

lik.vcv <- make.ou(phy, states, states.sd=SE,
                   control=list(method="vcv"))
fit.vcv <- find.mle(lik.vcv, c(1,0.1, mean(states)))
