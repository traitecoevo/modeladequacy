\name{var.contrast}
\alias{var.contrast}
\title{Variance of contrasts}
\usage{
  var.contrast(unit.tree)
}
\arguments{
  \item{unit.tree}{a \code{unit.tree} object}
}
\value{
  the estimate of the variance of the absolute value of the
  contrasts
}
\description{
  Calculates the variance of the absolute value of the
  contrasts
}
\details{
  This summary is used to evaluate whether the model is
  sufficiently capturing the variation in rates across the
  tree.

  The variance of contrasts is included as a default
  summary statistic in the function
  \code{\link{summ.stats}}
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## estimate variance of contrasts
var.contrast(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}}
}

