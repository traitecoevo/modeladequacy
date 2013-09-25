library(arbutus)
library(nlme)

## simulate data
set.seed(1)
t <- sim.bdtree(stop="taxa")
a <- sim.char(t, 1, model="BM")[,,]
b <- sim.char(t, 1, model="BM")[,,]

dat <- cbind.data.frame(a,b)

gls.bm <- gls(a~b, data=dat, correlation=corBrownian(phy=t))
gls.ou <- gls(a~b, data=dat, correlation=corMartins(value=1, phy=t))
gls.lam <- gls(a~b, data=dat, correlation=corPagel(value=1, phy=t))
gls.eb <- gls(a~b, data=dat, correlation=corBlomberg(value=0.1, phy=t, fixed=TRUE))

r.bm <- as.numeric(residuals(gls.bm))
names(r.bm) <- t$tip.label

r.ou <- as.numeric(residuals(gls.ou))
names(r.ou) <- t$tip.label

r.lam <- as.numeric(residuals(gls.lam))
names(r.lam) <- t$tip.label

r.eb <- as.numeric(residuals(gls.eb))
names(r.eb) <- t$tip.label

## fit using fit Continuous
fit.bm <- fitContinuous(phy=t, dat=r.bm, model="BM")
fit.ou <- fitContinuous(phy=t, dat=r.ou, model="OU")
fit.lam <- fitContinuous(phy=t, dat=r.lam, model="lambda")
fit.eb <- fitContinuous(phy=t, dat=r.eb, model="EB", bounds = list(a=c(mn=0.1, mx=0.11)))

## get parameter estimates and compare to estimates from gls
fit.bm$opt$sigsq
arbutus:::modelpars.gls(gls.bm)$sigsq

fit.ou$opt$alpha
arbutus:::modelpars.gls(gls.ou)$alpha
fit.ou$opt$sigsq
arbutus:::modelpars.gls(gls.ou)$sigsq

fit.lam$opt$lambda
arbutus:::modelpars.gls(gls.lam)$lambda
fit.lam$opt$sigsq
arbutus:::modelpars.gls(gls.lam)$sigsq

fit.eb$opt$sigsq
arbutus:::modelpars.gls(gls.eb)$sigsq

