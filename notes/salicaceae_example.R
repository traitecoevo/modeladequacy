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
lik <- make.ou(phy, states, states.sd=SE, control=list(method="pruning"))
dt <- find.mle(lik, x.init=c(1,0.1, mean(states)))

dt$par

## my result


#           s2         alpha         theta 
# 1.155839e-11  2.273490e+00 -8.623042e+03 

