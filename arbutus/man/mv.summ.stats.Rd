\name{mv.summ.stats}
\alias{mv.summ.stats}
\title{Multivariate measure of model adequacy}
\usage{
  mv.summ.stats(x)
}
\arguments{
  \item{x}{a \code{phy.ss} object inherited from
  \code{\link{compare.summ.stats}}}
}
\value{
  the Mahalanobis distance between the observed and
  simulated summary statistic.
}
\description{
  Computes Mahalanobis distance between the observed
  summary statistics and the simulated summary statistics
  as a multivariate measure of model fit
}
\details{
  This function computes the Mahalanobis distance between
  the observed and simulated summary statistics. The
  Mahalanobis distance (see
  \code{\link[stats]{mahalanobis}} is a unit-less and
  scale-invariant of the distance between a single data
  point (our observed summary statistics) and a common
  point (here the mean of the simulated summary
  statistics), taking into account the covariance between
  the summary statistics from the simulated data. It
  assumes that the distribution of summary statistics is
  multivariate normal. For the default summary statistics
  (see \code{\link{def.summ.stats}}, this condition should
  be met -- the exception being the
  \code{\link{ks.contrast}} statistic, which as it is
  bounded at 0, will not be. As a result, if
  \code{\link{ks.contrast}} is included in the set of
  summary statistic, this function will take the natural
  log of the values before computing the Mahalanobis
  distance. All other summary statistics will be assumed to
  be normally distributed and used as is.

  While the Mahalanobis distance may be a useful summary
  measure in some circumstances, we recommend checking the
  summary statistics individually for a number of reasons.
  First, our procedure for calculating the p-values for the
  summary statistics is general and does not depend on
  assumptions regarding the distribution of values. Second,
  the interpretation from the individual p-values is much
  more clear from the perspective of either posterior
  predictive or parametric bootstrapping theory. Third, and
  most importantly, the fact that some summary statistics
  capture the variation in the data better than others
  provides useful information as to how and why the model
  is inadequate.

  If only one set of observed summary statistics are
  available (e.g. from fitting a model using maximum
  likelihood to a single tree), a single distance will be
  returned. If multiple sets of summary statistics are
  available (e.g. from fitting a model using a Bayesian
  MCMC), the function will return a distribution of
  distances.
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## calculate default summary stats on observed data
ss.obs <- summ.stats(unit.tree, stats=NULL)

## simulate data on unit.tree
sims <- sim.char.unit(unit.tree, nsim=10)

## calculate default summary stats on simulated data
ss.sim <- summ.stats(sims, stats=NULL)

## compare simulated to observed summary statistics
res <- compare.summ.stats(ss.obs, ss.sim)

## calculate Mahalanobis distance
mv.summ.stats(res)
}
\seealso{
  \code{\link{compare.summ.stats}},
  \code{\link[stats]{mahalanobis}}
}

