\name{phy.model.check}
\alias{phy.model.check}
\title{Assess the adequacy of continuous trait models}
\usage{
  phy.model.check(x, nsim = 1000, stats = NULL, ...)
}
\arguments{
  \item{x}{a fitted model object or a
  \code{phylo}/\code{multiPhylo} (see Details).}

  \item{nsim}{the number of datasets to simulate. This is
  passed as an argument to the function
  \code{\link{sim.char.unit}}.}

  \item{stats}{a named list of summary statistics to
  calculate on observed and simulated datasets. See
  \code{\link{summ.stats}} for details.}

  \item{...}{additional arguments to be passed to
  \code{\link{as.unit.tree}}}
}
\description{
  Use summary statistics to assess model adequacy for
  phylogenetic models of continuous character evolution.
  This function is a simply a wrapper which does all the
  steps involved in evaluating model adequacy using the
  approach outlined in Pennell et al. (in prep.)
}
\details{
  \code{arbutus} was developed to evaluate the adequacy (or
  absolute goodness-of-fit) of phylogenetic models of
  continuous character evolution. The basic principle
  underlying the approach is that if the generating model
  is a Brownian motion process, the contrasts (sensu
  Felsenstein, 1985) will be independent and identically
  distributed (I.I.D.). We can evaluate this condition by
  calculating a set of summary statistics the contrasts on
  our observed data, then simulating datasets under a
  Brownian motion process along the phylogeny and
  calculating the summary statistics on the contrasts of
  each simulated data set. We can then compare our observed
  summary statistics to the simulated distribution of
  summary statistics.

  The I.I.D. property of the contrasts will hold only for
  the case of Brownian motion. However, the approach can be
  extended to any arbitrarily complex model of continuous
  character evolution (as long as the model is based on a
  multivariate normal distribution) by first rescaling the
  phylogeny based on the fitted parameter estimates,
  creating what we refer to as a 'unit.tree'. Once the
  rescaling is done, the contrasts will be I.I.D. if the
  model that was fit is a adequate one.

  P-values represent the two-tailed probability that the
  observed summary statistic (or distribution of summary
  statistics) came from the same distribution as the
  simulated summary statistics. Low p-values provide
  evidence that the model is inadequate.

  The function \code{phy.model.check} wraps several other
  functions: \enumerate{ \item{\code{\link{as.unit.tree}}:
  }{rescales phylogeny based on fitted parameter values}
  \item{\code{\link{summ.stats}}: }{calculates summary
  statistics on observed data}
  \item{\code{\link{sim.char.unit}}: }{simulates datasets
  under Brownian motion} \item{\code{\link{summ.stats}}:
  }{calculates summary statistics on simulated data}
  \item{\code{\link{compare.summ.stats}}: }{compares
  observed to simulated statistics} }

  \code{phy.model.check} interacts with objects produced by
  fitting evolutionary models using a variety of packages.
  Currently supported objects are as follows:

  \itemize{ \item a \code{gfit} object returned from
  fitting a model of continuous character evolution using
  \code{fitContinuous} in the \code{geiger} package.

  \item a \code{fit.mle} object returned from fitting a
  model of continuous character evolution using
  \code{find.mle} in the \code{diversitree} package. As the
  \code{fit.mle} object does not include all of the
  information required for creating a \code{unit.tree}, a
  second argument \code{lik} needs to be supplied,
  providing the likelihood function used in
  \code{find.mle}.

  \item a \code{mcmcsamples} object returned from fitting a
  model of continuous character evolution using MCMC
  methods in the \code{diversitree} package.
  \code{as.unit.tree} will apply the same trait dataset to
  a set of unit trees based on sampled parameters. By
  default this will create a unit tree for every sample in
  the mcmc chain. To modify this, additional arguments can
  be use.  \code{burnin} specifies how many samples to
  remove from the beginning of the chain.  \code{thin}
  specifies the thinning interval (e.g. if \code{thin=5},
  the function will create a unit tree from every fifth
  parameter set sampled.  \code{sample} specifies how many
  samples to draw from the MCMC run.  As the
  \code{mcmcsamples} object does not include all of the
  information required for creating a \code{unit.tree}, a
  second argument \code{lik} needs to be supplied,
  providing the likelihood function used in
  \code{mcmcsamples}.

  \item a \code{gls} object returned from fitting a
  phylogenetic generalized least squares model of character
  correlation using \code{gls} in the \code{nlme} package.

  \item a \code{pgls} object returned from fitting a
  phylogenetic generalized least squares model of character
  correlation using \code{pgls} in the \code{caper}
  package.

  \item a \code{phylolm} object returned from fitting a
  phylogenetic generalized linear model of character
  correlation using \code{phylolm} in the \code{phylolm}
  package. As the phylogeny is not returned with the
  \code{phylolm} object, a \code{phy} argument must also be
  specified.

  \item a \code{phylo} object. If a \code{phylo} object is
  supplied, the tree is assumed to have been rescaled
  previously. A second argument \code{data} must also be
  provided included the trait data as a named vector with
  names equal to the tip.labels of the phylogeny.

  \item a \code{multiPhylo object}. If a \code{multiPhylo}
  object is supplied, the tree is assumed to have been
  rescaled previously. A second argument \code{data} must
  also be provided included the trait data as a named
  vector with names equal to the tip.labels of the
  phylogenies. Note that this function will append the same
  data set to every tree in the \code{multiPhylo} object. }
}
\examples{
## finch data
data(finch)
phy <- finch$phy
data <- finch$data[,"wingL"]

\dontrun{
require(geiger)
## fit Brownian motion model
## using geiger's fitContinuous function
fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
                                 control=list(niter=10))

## check adequacy of BM model
## get p-values for default summary statistics
modelad.bm <- phy.model.check(fit.bm, nsim=10)
modelad.bm


## fit Ornstein-Uhlenbeck model
## again, using geiger's fitContinuous function
fit.ou <- suppressWarnings(fitContinuous(phy=phy, dat=data, model="OU",
                                                  control=list(niter=10)))

## check adequacy of OU model
modelad.ou <- phy.model.check(fit.ou, nsim=10)
}
}
\seealso{
  \code{\link{as.unit.tree}}, \code{\link{summ.stats}},
  \code{\link{sim.char.unit}},
  \code{\link{compare.summ.stats}}
}

