\name{plot.phy.ss}
\alias{plot.phy.ss}
\title{Plot summary statistics for model adequacy}
\usage{
  \method{plot}{phy.ss} (x, colour = NULL, ...)
}
\arguments{
  \item{x}{a 'phy.ss' object inherited from
  \code{\link{phy.model.check}} (or, alternatively
  \code{\link{compare.summ.stats}})}

  \item{colour}{a vector of two colours to be used in the
  plot; the first is for the observed summary statistics
  and the second for the simulated summary statistics.}

  \item{...}{additional arguments to be passed to
  \code{plot}. This is currently non-functional.}
}
\description{
  Plots the distribution of simulated summary statistics
  along with the observed summary statistics (either a
  point estimate or a distribution of values) for visual
  inspection of model adequacy
}
\details{
  This function is designed to give researchers a quick
  visual inspection of model adequacy by plotting the
  simulated distribution of summary statistics along with
  the empirical estimate(s). If only a point estimated of
  the observed summary statistics is input (such as from
  estimating the parameters on a single tree using maximum
  likelihood), the function will plot the simulated summary
  statistics with a vertical line representing the observed
  summary statistic. If a distribution of observed summary
  statistics is used, the function will plot both
  distributions (simulated and observed) on the same plot.
  This function is not designed to be flexible. If
  researchers want to create their own plots in a different
  style, all the necessary information is available in the
  'phy.ss' object output by both
  \code{\link{phy.model.check}} and
  \code{\link{compare.summ.stats}} and can be indexed by
  \code{x$summ.stats.obs} and \code{x$summ.stats.sim}.

  This plotting function using the
  \code{\link{ggplot2::ggplot2}} library.
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]

p <- phy.model.check(phy, data=dat, nsim=20)

plot(p)
}
\seealso{
  \code{\link{phy.model.check}},
  \code{\link{compare.summ.stats}}
}

