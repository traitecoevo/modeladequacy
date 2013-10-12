\name{sigsq.reml}
\alias{sigsq.reml}
\title{REML estimate of BM parameter}
\usage{
  sigsq.reml(unit.tree)
}
\arguments{
  \item{unit.tree}{a 'unit.tree' object}
}
\value{
  estimate of \eqn{sigma^2}
}
\description{
  Calculates the mean of the squared contrasts, which is
  equal to the REML estimate of the BM rate parameter. Used
  as a summary statistic to evaluate model adequacy
}
\details{
  This summary statistic is used to assess the ability of
  the model to adequate estimate the rate of evolution. On
  datasets simulated on a unit.tree (see
  \code{\link{sim.char.unit}}), the expected value is equal
  to 1.

  The REML estimate of \eqn{sigma^2} is included as a
  default summary statistic in the function
  \code{\link{summ.stats}}
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## Estimate sigsq from unit.tree

sigsq.reml(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}}
}

