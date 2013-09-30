\name{compare.summ.stats}
\alias{compare.summ.stats}
\title{Compare observed and simulated summary statistics}
\usage{
  compare.summ.stats(summ.stats.obs, summ.stats.sim)
}
\arguments{
  \item{summ.stats.obs}{a data.frame containing the
  observed summary statistics with the column names
  denoting the names of the statistics used.}

  \item{summ.stats.sim}{a data.frame containing the
  simulated summary statistics with the column names
  denoting the names of the statistics used.}
}
\value{
  an object of class 'phy.ss' which is a list containing
  the following objects: \itemize{ \item{\code{p.values}:
  }{Two--tailed p-values for all summary statistics}
  \item{\code{summ.stats.obs}: }{The input data frame of
  observed summary statistics} \item{\code{summ.stats.sim}:
  }{The input data frame of simulated summary statistics} }
}
\description{
  Calculates two--tailed p-values from comparing observed
  to simulated summary statistics in order to evaluate
  model adequacy for phylogenetic models of continuous
  character evolution.
}
\details{
  The column names must be the same for the observed and
  simulated summary statistics (though not necessarily in
  the same order). Both \code{summ.stats.obs} and
  \code{summ.stats.sim} are designed to be obtained from
  the function \code{\link{summ.stats}}.

  If there is only 1 row in \code{summ.stats.obs} (such as
  the case if ML estimate of model is used), then the
  function will compare the observed summary statistic(s)
  to the distribution of simulated simulated summary
  statistics.

  If there is more than 1 row in \code{summ.stats.obs}
  (such as the case if we have a posterior distribution of
  parameter estimates or the summary statistics are
  calculated across many trees), then the function requires
  that the number of rows in \code{summ.stats.obs} and
  \code{summ.stats.sim} be equal. This is because the
  function calculates the p-values by performing pairwise
  comparisons. the function will automatically check that
  this condition is met and produce an error if it is not.

  In both cases, the p-values that are returned are
  two--tailed p-values as the function considers the
  probability of the observed statistics being larger or
  smaller than the simulated statistics.

  The distribution of simulated summary statistics along
  with the observed summary statistic using a generic
  \code{\link{plot}} function (see
  \code{\link{plot.phy.ss}}) for more details. For the
  purposes of plotting, the observed and simulated summary
  statistics used as arguments are returned with the
  p-values.

  The p-values can be extracted from the output of
  \code{compare.summ.stats} using the function
  \code{\link{pval.summ.stats}}
}
\examples{
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
dat <- td$data[,"wingL"]
unit.tree <- as.unit.tree(phy, dat)

## calculate default summary stats on observed data
ss.obs <- summ.stats(unit.tree, stats=NULL)

## simulate data on unit.tree
sims <- sim.char.unit(unit.tree, nsim=10)

## calculate default summary stats on simulated data
ss.sim <- summ.stats(sims, stats=NULL)

## compare simulated to observed summary statistics
res <- compare.summ.stats(ss.obs, ss.sim)

## extract p-values
pval.summ.stats(res)

## plot results
plot(res)
}
\author{
  Matt Pennell
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{plot.phy.ss}},
  \code{\link{pval.summ.stats}},
  \code{\link{phy.model.check}}
}

