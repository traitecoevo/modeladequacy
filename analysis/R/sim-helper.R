library(diversitree)
library(arbutus)
library(multicore)

## use interior functions from arbutus
sim_bm <- arbutus:::sim.char.std.bm
model_phylo_bm <- arbutus:::model_phylo_bm
model_phylo_ou <- arbutus:::model_phylo_ou
model_phylo_eb <- arbutus:::model_phylo_eb

## set path
path.sim <- function()
    "output/sims/"

## overwrite function from model-adequacy-fits to change output
## Now writes to file

model_ad_simfit <- function(phy, states, pars, model, type, seed=1) {

  model <- match.arg(model, c("BM", "OU", "EB"))
  type  <- match.arg(type,  c("ml", "bayes"))
  SE    <- pars$SE 

  # Make the analyses recomputable by using the same seed each time:
  set.seed(seed)

  make.lik <- switch(model, BM=make.bm, OU=make.ou, EB=make.eb)
  lik <- make.lik(phy, states, SE,
                  control=list(method="pruning", backend="C"))

  # Start at REML estimate of sigsq for all models, using arbutus'
  # function to estiate this.
  s2 <- pic_stat_msig(make_unit_tree(phy, data=drop(states)))

  # ML bounds, starting points and priors (priors only used in the
  # Bayesian analysis)
  lower <- 0
  upper <- Inf
  if (model == "BM") {
    start <- s2
    prior <- make.prior.bm(s2.lower=0, s2.upper=2)
  } else if (model == "OU") {
    upper <- bounds.ou(phy, states)
    start <- c(s2, 0.05)
    prior <- make.prior.ou(s2.lower=0, s2.upper=2,
                           ln.mean=log(0.5), ln.sd=log(1.5))
  } else if (model == "EB") {
    lower <- c(0, -1)
    start <- c(s2, -0.1)
    prior <- make.prior.eb(s2.lower=0, s2.upper=2,
                           a.lower=-1, a.upper=0)
  }

  if (type == "ml") {
    fit <- find.mle(lik, x.init=start)
    ic  <- AIC(fit)
    ic.name <- "aic"

    # Assess adequacy of all models
    ma <- arbutus(fit)
  } else if (type == "bayes") {
    # Some general parameters
    pilot  <- 100
    burnin <- 1000
    nsteps <- 10000

    # Run short chain to obtain appropriate step size for MCMC
    tmp <- mcmc(lik, x.init=start, nsteps=pilot, prior=prior, w=1,
                print.every=0)

    w <- diff(apply(coef(tmp), 2, range))

    # Full chain:
    samples <- mcmc(lik, x.init=start, nsteps=nsteps, prior=prior, w=w,
                    print.every=0)

    ic <- dic.mcmcsamples(samples, burnin=burnin)
    ic.name <- "dic"

    # Assess adequacy of all models    
    ma <- arbutus(samples, burnin=burnin, sample=1000)
  }

  p <- as.numeric(ma$p.values)
  res <- list()
  res$model <- model
  res$size <- Ntip(phy)
  res$sigsq <- pars$sigsq
  res$alpha <- pars$alpha
  res$a <- pars$a
  res$SE <- SE
  res$m.sig <- p[1]
  res$c.var <- p[2]
  res$s.var <- p[3]
  res$s.asr <- p[4]
  res$s.hgt <- p[5]
  res$d.cdf <- p[6]
  res
}


## Functions for simulation and fitting data
model_ad_sim <- function(pars, model, type, n.taxa){
    phy <- tree.bd(c(1,0), max.taxa=n.taxa)
    ## rescale tree to be unit length
    phy$edge.length <- phy$edge.length / max(branching.times(phy))

    ## simulate data by constructing unit tree
    ## ignore SE for time being
    tmp.p <- pars
    tmp.p$SE <- 0

    states <- switch(model,
                     BM=sim_bm(model_phylo_bm(phy, tmp.p)),
                     OU=sim_bm(model_phylo_ou(phy, tmp.p)),
                     EB=sim_bm(model_phylo_eb(phy, tmp.p)))
        
    ## add error
    ## NOTE: Need to check SE vs. SD!!
    SE <- pars$SE
    states <- states[,1] + rnorm(n.taxa, sd=SE)
    model_ad_simfit(phy, states, pars, model, type) 
}


## Wrapper function to simulate multiple times
model_ad_sim_multi <- function(pars, model, type, n.taxa, n.sims, write=TRUE, filename){
    res <- mclapply(seq_len(nsims), function(x)
                    model_ad_sim(pars, model, type, n.taxa),
                    mc.cores=mc.cores, mc.preschedule=FALSE)
    res <- do.call(rbind, res)

    if (write)
        saveRDS(res, paste(path.sim(), filename, ".rds", sep=""))

    res 
}



