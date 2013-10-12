\name{ks.pic}
\alias{ks.pic}
\title{Evaluate normality of contrasts}
\usage{
  ks.pic(unit.tree)
}
\arguments{
  \item{unit.tree}{a unit.tree object}
}
\value{
  the D-statistic from a KS-test
}
\description{
  Performs Kolmogorov-Smirnov test comparing distribution
  of contrasts with that of a normal distribution with mean
  0 and variance equal to the square root of the mean of
  squared contrasts
}
\details{
  This summary statistics is used to evaluate whether the
  assumption of multivariate normailty is appropriate. If
  the model which generated the data is the fitted model,
  we expect the square root ofthe mean of squared contrasts
  to be equal to 1. The empirical estimate is used rather
  than assume a variance of 1 to reduce the overlap between
  this summary statistic and the REML estimate of
  \eqn{sigma^2} (see \code{\link{reml.sigsq}}). The
  Kolmogorov-Smirnov (KS) test is a non-parameteric test
  which computes the maximum distance \eqn{D} between two
  cumulative distribution functions. Running the test
  multiple times on the same data will produce slightly
  different values due to the fact that the null
  distribution is produced by randomly drawing from a
  normal distribution.

  The KS-D statistic is included as a default summary
  statistic in the function \code{\link{summ.stats}}.
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## KS-D statistic
ks.pic(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}},
  \code{\link{stats::ks.test}}
}

