## This is needed until I nail down some issues with argument passing
## and peering into lists:
run_model_ad_ml <- function(dat) {
  lapply(dat, run_model_ad, "ml")
  ## mclapply(dat, run_model_ad, "ml", mc.preschedule=FALSE)
}
run_model_ad_bayes <- function(dat) {
  mclapply(dat, run_model_ad, "bayes", mc.preschedule=FALSE)
}

run_model_ad_ml_noerr <- function(dat) {
  mclapply(dat, run_model_ad, "ml", err=FALSE, mc.preschedule=FALSE)
}
run_model_ad_ml_nolog <- function(dat) {
  mclapply(dat, run_model_ad, "ml", log=FALSE, mc.preschedule=FALSE)
}


# Little wrapper function that loads data, fits a model (ML or MCMC),
# checks the model adequacy, and then saves the results in a file.  If
# the output file exists, this is skipped.
run_model_ad <- function(dat, type, verbose=TRUE, err=TRUE, log=TRUE) {
  if (verbose) {
    message(sprintf("%s: %s / %s / %s",
                    type, dat$trait, dat$rank, dat$taxa))
  }
  models <- c("BM", "OU", "EB")
  res <- lapply(models, function(m) model_ad(dat, m, type, err=err, log=log))
  names(res) <- models
  res
}

# This function does all the ML or Bayesian analysis, depending on the
# 'type' argument.  It takes a data set (which contains a tree,
# character states, and standard errors for the character states),
# fits a model using ML and then assesses the adequacy of that model.
#
# Most of the ugliness is putting (hopefully) sensible bounds on the
# parameters to avoid ridges in likelihood space that can happen with
# OU and EB.
#
# Note that the entire adequacy part of the function is the line
#      ma <- arbutus(fit)
# or
#     ma <- arbutus(samples)
model_ad <- function(data, model, type, seed=1, err=TRUE, log=TRUE) {
  model <- match.arg(model, c("BM", "OU", "EB"))
  type  <- match.arg(type,  c("ml", "bayes"))
  ## Extract components from the pre-prepared data object.
  phy    <- data$phy

  ## Different options for the analysis:
  if (log) {
    states <- drop(data$states_log)
    se     <- data$se_log
  } else {
    states <- drop(data$states_raw)
    se     <- data$se_raw
  }
  if (!err) {
    se <- 0.0
  }
  
  # Make the analyses recomputable by using the same seed each time:
  set.seed(seed)

  make_lik <- switch(model, BM=make.bm, OU=make.ou, EB=make.eb)
  lik <- make_lik(phy, states, se,
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
    prior <- make_prior_bm(s2_lower=0, s2_upper=2)
  } else if (model == "OU") {
    upper <- bounds_ou(phy, states)
    start <- c(s2, 0.05)
    prior <- make_prior_ou(s2_lower=0, s2_upper=2,
                           ln_mean=log(0.5), ln_sd=log(1.5))
  } else if (model == "EB") {
    lower <- c(0, -1)
    start <- c(s2, -0.1)
    prior <- make_prior_eb(s2_lower=0, s2_upper=2,
                           a_lower=-1, a_upper=0)
  }

  if (type == "ml") {
    fit <- find.mle(lik, x.init=start)
    ic  <- AIC(fit)
    ic_name <- "aic"

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

    ic <- dic_mcmcsamples(samples, burnin=burnin)
    ic_name <- "dic"

    # Assess adequacy of all models    
    ma <- arbutus(samples, burnin=burnin, sample=1000)
  }

  # For all cases, the general information about three we want is the same:
  info       <- data[c("taxa", "rank", "trait")]
  info$size  <- Ntip(phy)
  info$age   <- max(branching.times(phy))
  info$model <- model
  info$type  <- type
  info[[ic_name]] <- ic

  # Returning information about the tree, the assessment of model
  # adequacy and the fit (either the ML object or the MCMC samples).
  list(info=info, ma=ma, fit=if(type == "ml") fit else samples)
}

## function for making prior for bm
## just a wrapper for make.prior.uniform
make_prior_bm <- function(s2_lower, s2_upper) {
  make.prior.uniform(lower=s2_lower, upper=s2_upper)
}

## define function to create a lognormal prior
## to be used in make_prior_ou
make_prior_lognormal <- function(ln_mean, ln_sd) {
  function(x) sum(dlnorm(x, meanlog=ln_mean, sdlog=ln_sd, log=TRUE))
}

## function for making ou prior
## s2 is given a uniform prior
## alpha is given a lognormal prior
make_prior_ou <- function(s2_lower, s2_upper, ln_mean, ln_sd){
  p_s2 <- make.prior.uniform(lower=s2_lower, upper=s2_upper)
  p_al <- make_prior_lognormal(ln_mean=ln_mean, ln_sd=ln_sd)

  function(pars) {
    p_s2(pars[[1]]) + p_al(pars[[2]])
  }
}

## function for making eb prior
## s2 is given a uniform prior
## a is givne a uniform prior
## the prior on s2 should be different than that of a
## s2 > 0 and a < 0
make_prior_eb <- function(s2_lower, s2_upper, a_lower, a_upper){
  p_s2 <- make.prior.uniform(lower=s2_lower, upper=s2_upper)
  p_a  <- make.prior.uniform(lower=a_lower, upper=a_upper)

  function(pars) {
    p_s2(pars[[1]]) + p_a(pars[[2]])
  }
}

## create function for setting OU bounds for ML analyses
## this is the only one where bounds to be necessary
bounds_ou <- function(phy, states){
  ## find the shortest terminal branch
  sht <- min(phy$edge.length[phy$edge[,2] <= Ntip(phy)])
  ## upper bound for alpha: greater than shortest branch
  ## upper bound for for sigsq: set to 2
  c(2, log(2)/sht)
}
